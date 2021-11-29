"use strict";


var viewport = {
    iRz: 0.9,
    iRx: 0.2,
    iSc: 1.3,
    renderNeeded: true
};

var texture = {
    id: -1,
    object: {},
    texture: null,
    dims: [1, 1, 1],
    bbox: [0, 0, 0]
};

var requestCache = {};

function loadShaderSource(path) {
    var source = "";
    if (requestCache[path] != undefined) {
        source = requestCache[path];
    }
    else {
        var request = new XMLHttpRequest();
        request.open("GET", path, false);
        request.send(null);
        if (request.status == 200) {
            source = request.responseText;
            requestCache[path] = source;
        }
    }
    return source;
}

function loadTexture(gl, url) {

    function onload(treeBuffer) {
        const tex = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, tex);

        console.log(treeBuffer.length);
        const BUFFER_SIZE = 1024;
        var treeBufferPadded = new Uint8Array(BUFFER_SIZE * BUFFER_SIZE * 4);
        treeBufferPadded.set(treeBuffer, 0);

        console.log(treeBufferPadded.length);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8UI,
            1024, 1024, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_BYTE,
            treeBufferPadded);

        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_BASE_LEVEL, 0);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAX_LEVEL, 0);

        texture.texture = tex;

        viewport.renderNeeded = true;
    }

    function onerror() {
        document.getElementById("volume-dims").innerHTML
            = "<span style='color:red;'>Failed to load volume.</span>";
        document.getElementById("volume-select").disabled = false;
    }

    if (requestCache[url] != undefined) {
        onload(requestCache[url]);
        return;
    }

    var req = new XMLHttpRequest();
    req.open("GET", url, true);
    req.responseType = "arraybuffer";

    req.onload = function (e) {
        if (req.status == 200) {
            var treeBuffer = new Uint8Array(req.response);
            requestCache[url] = treeBuffer;
            onload(treeBuffer);
        }
        else {
            onerror();
        }
    };
    req.onerror = function (e) {
        onerror();
    };
    req.send();
}


// initialize WebGL: load and compile shader, initialize buffers
function initWebGL(gl) {

    console.time("request glsl code");
    var vsSource = `#version 300 es
        precision highp float;

        in vec4 vertexPosition;

        void main(void) {
            gl_Position = vertexPosition;
            return;
        }
    `;
    var fsSource = loadShaderSource("frag.glsl");
    console.timeEnd("request glsl code");

    // compile shaders
    function initShaderProgram(gl, vsSource, fsSource) {
        function loadShader(gl, type, source) {
            var shader = gl.createShader(type); // create a new shader
            gl.shaderSource(shader, source); // send the source code to the shader
            gl.compileShader(shader); // compile shader
            if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) // check if compiled succeed
                throw new Error(gl.getShaderInfoLog(shader)); // compile error message
            return shader;
        }
        var vShader = loadShader(gl, gl.VERTEX_SHADER, vsSource);
        var fShader = loadShader(gl, gl.FRAGMENT_SHADER, fsSource);
        // create the shader program
        var shaderProgram = gl.createProgram();
        gl.attachShader(shaderProgram, vShader);
        gl.attachShader(shaderProgram, fShader);
        gl.linkProgram(shaderProgram);
        // if creating shader program failed
        if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS))
            throw new Error(gl.getProgramInfoLog(shaderProgram));
        return shaderProgram;
    }
    console.time("compile shader");
    var shaderProgram = initShaderProgram(gl, vsSource, fsSource);
    console.timeEnd("compile shader");

    // position buffer
    var positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    var positions = [-1, 1, 1, 1, -1, -1, 1, -1];
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(positions), gl.STATIC_DRAW);

    // return a JSON object
    var programInfo = {
        program: shaderProgram,
        attribLocations: { // attribute variables, receive values from buffers
            vertexPosition: gl.getAttribLocation(shaderProgram, 'vertexPosition'),
        },
        uniformLocations: { // uniform variables
            iRz: gl.getUniformLocation(shaderProgram, 'iRz'),
            iRx: gl.getUniformLocation(shaderProgram, 'iRx'),
            iSc: gl.getUniformLocation(shaderProgram, 'iSc'),
            iResolution: gl.getUniformLocation(shaderProgram, "iResolution"),
            uTreeBuffer: gl.getUniformLocation(shaderProgram, "uTreeBuffer"),
        },
        buffers: {
            positionBuffer: positionBuffer,
        },
    };
    return programInfo;
}


// call this function to re-render
function drawScene(gl, programInfo) {

    // clear the canvas
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clearColor(0, 0, 0, 1);
    gl.clearDepth(-1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.GEQUAL);

    // tell WebGL how to pull out the positions from the position buffer into the vertexPosition attribute
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, programInfo.buffers.positionBuffer);
        gl.vertexAttribPointer(
            programInfo.attribLocations.vertexPosition,
            numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(programInfo.attribLocations.vertexPosition);
    }

    // indice buffer
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, programInfo.buffers.indiceBuffer);

    // make sure it uses the program
    gl.useProgram(programInfo.program);

    // set shader uniforms
    // https://webglfundamentals.org/webgl/lessons/webgl-shaders-and-glsl.html
    gl.uniform1f(programInfo.uniformLocations.iRz, viewport.iRz + 1e-4);
    gl.uniform1f(programInfo.uniformLocations.iRx, viewport.iRx + 1e-4);
    gl.uniform1f(programInfo.uniformLocations.iSc, viewport.iSc);
    gl.uniform2f(programInfo.uniformLocations.iResolution, canvas.clientWidth, canvas.clientHeight);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, texture.texture);
    gl.uniform1i(programInfo.uniformLocations.uTreeBuffer, 0);

    // render
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }
}


function main() {
    const canvas = document.getElementById("canvas");
    const gl = canvas.getContext("webgl2") || canvas.getContext("experimental-webgl2");
    if (gl == null) throw ("Error: `canvas.getContext(\"webgl2\")` returns null. Your browser may not support WebGL 2.");

    // load WebGL
    var programInfo = initWebGL(gl);
    loadTexture(gl, "tree1.bin");

    // rendering
    let then = 0;
    function render_main(now) {
        if (viewport.renderNeeded) {
            // display fps
            now *= 0.001;
            var time_delta = now - then;
            then = now;
            if (time_delta != 0) {
                document.getElementById("fps").textContent = (1.0 / time_delta).toFixed(1) + " fps";
            }

            canvas.width = canvas.style.width = window.innerWidth;
            canvas.height = canvas.style.height = window.innerHeight;
            drawScene(gl, programInfo);

            viewport.renderNeeded = false;
        }
        requestAnimationFrame(render_main);
    }
    requestAnimationFrame(render_main);

    // interactions
    canvas.addEventListener("wheel", function (e) {
        e.preventDefault();
        var sc = Math.exp(-0.0005 * e.wheelDeltaY);
        viewport.iSc *= sc;
        viewport.renderNeeded = true;
    }, { passive: false });
    var mouseDown = false;
    canvas.addEventListener("pointerdown", function (event) {
        //event.preventDefault();
        mouseDown = true;
    });
    window.addEventListener("pointerup", function (event) {
        event.preventDefault();
        mouseDown = false;
    });
    window.addEventListener("resize", function (event) {
        canvas.width = canvas.style.width = window.innerWidth;
        canvas.height = canvas.style.height = window.innerHeight;
        viewport.renderNeeded = true;
    });
    canvas.addEventListener("pointermove", function (e) {
        if (mouseDown) {
            viewport.iRx += 0.01 * e.movementY;
            viewport.iRz -= 0.01 * e.movementX;
            viewport.renderNeeded = true;
        }
    });
}

window.onload = function (event) {
    setTimeout(function () {
        try {
            main();
        } catch (e) {
            console.error(e);
            document.body.innerHTML = "<h1 style='color:red;'>" + e + "</h1>";
        }
    }, 0);
};