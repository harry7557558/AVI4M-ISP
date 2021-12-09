"use strict";


// request shader source
function loadShaderSource(path) {
    var source = "";
    var request = new XMLHttpRequest();
    request.responseType = "";
    request.open("GET", path, false);
    request.send(null);
    if (request.status == 200) {
        source = request.responseText;
    }
    else {
        throw ("Error loading shader source (" + request.status + "): " + path);
    }
    const include_regex = /\n\#include\s+[\<\"](.*?)[\>\"]/;
    while (include_regex.test(source)) {
        var match = source.match(include_regex);
        var include_source = loadShaderSource(match[1]);
        source = source.replace(match[0], "\n" + include_source + "\n");
    }
    return source;
}

// compile shader
function initShaderProgram(gl, vsSource, fsSource) {
    function loadShader(gl, type, source) {
        var shader = gl.createShader(type); // create a new shader
        gl.shaderSource(shader, source); // send the source code to the shader
        gl.compileShader(shader); // compile shader
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) // check if compiled succeed
            throw "Shader compile error:<br/>" + gl.getShaderInfoLog(shader); // compile error message
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

// load tree data structure for an object
function loadTreeTexture(renderer, url, texture_name) {
    const gl = renderer.gl;

    function onload(treeBuffer) {
        const tex = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, tex);

        console.log(treeBuffer.length);
        const BUFFER_SIZE = 1024;
        var treeBufferPadded = new Uint8Array(BUFFER_SIZE * BUFFER_SIZE * 4);
        treeBufferPadded.set(treeBuffer, 0);

        console.log(treeBufferPadded.length);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8UI,
            BUFFER_SIZE, BUFFER_SIZE, 0,
            gl.RGBA_INTEGER, gl.UNSIGNED_BYTE,
            treeBufferPadded);

        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_BASE_LEVEL, 0);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAX_LEVEL, 0);

        renderer.textures[texture_name] = tex;
        renderer.uniforms.iFrame = 0;
        renderer.renderNeeded = true;
    }

    var req = new XMLHttpRequest();
    req.open("GET", url, true);
    req.responseType = "arraybuffer";
    req.onerror = function (e) {
        alert("Failed to load texture " + texture_name);
    };
    req.onload = function (e) {
        if (req.status == 200) {
            var treeBuffer = new Uint8Array(req.response);
            onload(treeBuffer);
        }
        else {
            req.onerror();
        }
    };
    req.send();
}


// call this function to re-render
function drawScene(renderer) {
    let gl = renderer.gl;
    let programs = renderer.programs;
    let uniforms = renderer.uniforms;
    let buffers = renderer.buffers;
    let textures = renderer.textures;

    // clear the canvas
    gl.viewport(0, 0, uniforms.iResolution.x, uniforms.iResolution.y);
    gl.clearColor(0, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);

    // tell WebGL how to pull out the positions from the position buffer into the vertexPosition attribute
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, buffers.positionBuffer);
        var vertexPosition = gl.getAttribLocation(programs.renderProgram, 'vertexPosition');
        gl.vertexAttribPointer(
            vertexPosition,
            numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(vertexPosition);
    }

    // use the renderer program
    gl.useProgram(programs.renderProgram);

    // set shader uniforms
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iRz"), uniforms.iRz);
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iRx"), uniforms.iRx);
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iSc"), uniforms.iSc);
    gl.uniform2f(
        gl.getUniformLocation(programs.renderProgram, "iResolution"),
        uniforms.iResolution.x, uniforms.iResolution.y);
    gl.uniform1i(
        gl.getUniformLocation(programs.renderProgram, "iFrame"),
        uniforms.iFrame++);

    // set texture
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, textures.treeSampler);
    gl.uniform1i(gl.getUniformLocation(programs.renderProgram, "uTreeBuffer"), 0);

    // render
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }
}


function main() {

    // get WebGL context
    const canvas = document.getElementById("canvas");
    const gl = canvas.getContext("webgl2") || canvas.getContext("experimental-webgl2");
    if (gl == null) throw ("Failed to load WebGL2");

    var renderer = {
        canvas: canvas,
        gl: gl,
        extentions: {},
        programs: {
            renderProgram: null
        },
        uniforms: {
            iRz: 0.9 + 1e-4,
            iRx: 0.2 + 1e-4,
            iSc: 0.25,
            iResolution: { x: 1, y: 1 },
            iFrame: 0
        },
        buffers: {
            positionBuffer: null
        },
        textures: {
            treeSampler: null
        },
        statistics: {
            recordedTime: []
        },
        renderNeeded: true
    };
    renderer.extentions["EXT_disjoint_timer_query_webgl2"] = gl.getExtension("EXT_disjoint_timer_query_webgl2");
    renderer.uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
    renderer.uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;

    // load object - do this early because it takes time
    loadTreeTexture(renderer, "tree6.bin", "treeSampler");

    // load shaders
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
    console.time("compile shader");
    renderer.programs.renderProgram = initShaderProgram(gl, vsSource, fsSource);
    console.timeEnd("compile shader");

    // position buffer
    var positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    var positions = [-1, 1, 1, 1, -1, -1, 1, -1];
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(positions), gl.STATIC_DRAW);
    renderer.buffers.positionBuffer = positionBuffer;

    // rendering
    let then = 0;
    function render_main(now) {
        if (renderer.renderNeeded) {
            // display fps (time moving average)
            let recordedTime = renderer.statistics.recordedTime;
            var time_elapsed = 0.001 * (now - then);
            recordedTime.push(time_elapsed);
            renderer.statistics.recordedTime = recordedTime.slice(Math.max(recordedTime.length - 120, 0));
            var frame = 0, time = 0.0, maxTime = 1.0;
            for (var i = recordedTime.length - 1; i >= 0 && time < maxTime; i--) {
                frame += 1;
                time += recordedTime[i];
            }
            var avrfps = frame / time;
            document.getElementById("fps").innerHTML =
                renderer.uniforms.iResolution.x + " × " + renderer.uniforms.iResolution.y + "<br/>" +
                (1000.0 / avrfps).toFixed(0) + " ms, " + avrfps.toFixed(1) + " fps";
            drawScene(renderer);
            //renderer.renderNeeded = false;
        }
        then = now;
        if (renderer.textures.treeSampler == null) {
            document.getElementById("fps").innerHTML = "Loading model..."
        }
        requestAnimationFrame(render_main);
    }
    requestAnimationFrame(render_main);

    // interactions
    var mouseDown = false;
    canvas.addEventListener("wheel", function (e) {
        e.preventDefault();
        var sc = Math.exp(-0.0005 * e.wheelDeltaY);
        if (true) {
            renderer.uniforms.iSc *= sc;
            renderer.renderNeeded = true;
        }
    }, { passive: false });
    canvas.addEventListener("pointerdown", function (event) {
        //event.preventDefault();
        mouseDown = true;
    });
    window.addEventListener("pointerup", function (event) {
        event.preventDefault();
        mouseDown = false;
    });
    window.addEventListener("resize", function (event) {
        renderer.uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
        renderer.uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;
        renderer.renderNeeded = true;
    });
    canvas.addEventListener("pointermove", function (e) {
        if (mouseDown) {
            renderer.uniforms.iRx += 0.01 * e.movementY;
            renderer.uniforms.iRz -= 0.01 * e.movementX;
            renderer.renderNeeded = true;
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