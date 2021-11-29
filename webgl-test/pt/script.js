"use strict";


var uniforms = {
    iResolution: { x: 1, y: 1 },
    iFrame: 0,
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


function createSampleTexture(gl, width, height) {
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    const level = 0;
    const internalFormat = gl.RGBA32F;
    const border = 0;
    const format = gl.RGBA;
    const type = gl.FLOAT;
    const data = null;
    gl.texImage2D(gl.TEXTURE_2D, level, internalFormat,
        width, height, border,
        format, type, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    return tex;
}

function createRenderTarget(gl, width, height) {
    const tex = createSampleTexture(gl, width, height);
    const framebuffer = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    return {
        texture: tex,
        framebuffer: framebuffer
    };
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
    var fsSourceRender = loadShaderSource("frag.glsl");
    var fsSourceDisplay = `#version 300 es
        precision highp float;

        out vec4 fragColor;

        uniform sampler2D iChannel0;

        void main(void) {
            vec3 col = texelFetch(iChannel0, ivec2(gl_FragCoord.xy), 0).xyz;
            fragColor = vec4(col, 1.0);
        }
    `;
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
    console.time("compile shaders");
    var renderProgram = initShaderProgram(gl, vsSource, fsSourceRender);
    var displayProgram = initShaderProgram(gl, vsSource, fsSourceDisplay);
    console.timeEnd("compile shaders");

    // position buffer
    var positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    var positions = [-1, 1, 1, 1, -1, -1, 1, -1];
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(positions), gl.STATIC_DRAW);

    // frame buffer
    var renderTexture = createRenderTarget(gl, uniforms.iResolution.x, uniforms.iResolution.y);

    // sample texture
    var sampleTexture = createSampleTexture(gl, uniforms.iResolution.x, uniforms.iResolution.y);

    // return a JSON object
    var programInfo = {
        programs: {
            renderProgram: renderProgram,
            displayProgram: displayProgram
        },
        buffers: {
            positionBuffer: positionBuffer
        },
        textures: {
            renderTexture: renderTexture,  // texture + framebuffer
            sampleTexture: sampleTexture  // texture
        }
    };
    return programInfo;
}


// call this function to re-render
function drawScene(gl, programInfo) {

    /* RENDER TO FRAMEBUFFER */
    gl.bindFramebuffer(gl.FRAMEBUFFER, programInfo.textures.renderTexture.framebuffer);
    //gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, uniforms.iResolution.x, uniforms.iResolution.y);
    if (uniforms.iFrame == 0) {
        gl.clearColor(0, 1, 0.5, 0);
        gl.clearDepth(-1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }
    gl.useProgram(programInfo.programs.renderProgram);

    // set position buffer
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, programInfo.buffers.positionBuffer);
        var vertexPosition = gl.getAttribLocation(programInfo.programs.renderProgram, 'vertexPosition');
        gl.vertexAttribPointer(vertexPosition, numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(vertexPosition);
    }

    // set shader uniforms
    gl.uniform2f(
        gl.getUniformLocation(programInfo.programs.renderProgram, "iResolution"),
        uniforms.iResolution.x, uniforms.iResolution.y);
    gl.uniform1i(
        gl.getUniformLocation(programInfo.programs.renderProgram, "iFrame"),
        uniforms.iFrame++);

    // set texture
    gl.bindTexture(gl.TEXTURE_2D, programInfo.textures.sampleTexture);
    //gl.bindTexture(gl.TEXTURE_2D, programInfo.textures.renderTexture.texture);
    gl.uniform1i(gl.getUniformLocation(programInfo.programs.renderProgram, "iChannel0"), 0);

    // render to framebuffer
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }

    // copy rendered image
    gl.bindTexture(gl.TEXTURE_2D, programInfo.textures.sampleTexture);
    gl.copyTexImage2D(gl.TEXTURE_2D,
        0, gl.RGBA32F, 0, 0, uniforms.iResolution.x, uniforms.iResolution.y, 0);

    /* RENDERING TO DISPLAY */
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, uniforms.iResolution.x, uniforms.iResolution.y);
    gl.clearColor(0, 0.5, 1, 1);
    gl.clearDepth(-1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.useProgram(programInfo.programs.displayProgram);

    // set position buffer
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, programInfo.buffers.positionBuffer);
        var vertexPosition = gl.getAttribLocation(programInfo.programs.displayProgram, 'vertexPosition');
        gl.vertexAttribPointer(vertexPosition, numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(vertexPosition);
    }

    // set texture
    gl.bindTexture(gl.TEXTURE_2D, programInfo.textures.sampleTexture);
    //gl.bindTexture(gl.TEXTURE_2D, programInfo.textures.renderTexture.texture);
    gl.uniform1i(gl.getUniformLocation(programInfo.programs.displayProgram, "iChannel0"), 0);

    // render to canvas
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }
}


function main() {
    const canvas = document.getElementById("canvas");
    uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
    uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;

    const gl = canvas.getContext("webgl2") || canvas.getContext("experimental-webgl2");
    if (gl == null) throw ("Error: `canvas.getContext(\"webgl2\")` returns null. Your browser may not support WebGL 2.");
    var ext = gl.getExtension("EXT_color_buffer_float");
    if (!ext) throw ("Error: `gl.getExtension(\"EXT_color_buffer_float\")` returns null.");

    // load WebGL
    var programInfo = initWebGL(gl);

    // rendering
    let then = 0;
    function render_main(now) {
        // display fps
        now *= 0.001;
        var time_delta = now - then;
        then = now;
        if (time_delta != 0) {
            document.getElementById("fps").innerHTML =
                "Frame " + (uniforms.iFrame + 1) + "<br/>" +
                (1.0 / time_delta).toFixed(1) + " fps";
        }

        // render
        drawScene(gl, programInfo);

        requestAnimationFrame(render_main);
    }
    requestAnimationFrame(render_main);

    // interactions
    canvas.addEventListener("wheel", function (e) {
        e.preventDefault();
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
        // update resolution
        uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
        uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;
        // update framebuffer
        gl.deleteFramebuffer(programInfo.textures.renderTexture.framebuffer);
        gl.deleteTexture(programInfo.textures.renderTexture.texture);
        programInfo.textures.renderTexture = createRenderTarget(gl, uniforms.iResolution.x, uniforms.iResolution.y);
        // update copy buffer
        gl.deleteTexture(programInfo.textures.sampleTexture);
        programInfo.textures.sampleTexture = createSampleTexture(gl, uniforms.iResolution.x, uniforms.iResolution.y);
        uniforms.iFrame = 0;
    });
    canvas.addEventListener("pointermove", function (e) {
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