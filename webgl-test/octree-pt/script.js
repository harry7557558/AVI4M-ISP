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
    return source;
}

// compile shader
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

// load tree data structure for an object
function loadTreeTexture(renderer, url, texture_name) {
    const gl = renderer.gl;

    function onload(treeBuffer) {
        const tex = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, tex);

        const BUFFER_SIZE = 1024;
        var treeBufferPadded = new Uint8Array(BUFFER_SIZE * BUFFER_SIZE * 4);
        treeBufferPadded.set(treeBuffer, 0);
        console.log(treeBuffer.length, 4 * BUFFER_SIZE * BUFFER_SIZE);

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

// create render targets
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
    const sampler = createSampleTexture(gl, width, height);
    return {
        framebuffer: framebuffer,
        texture: tex,
        sampler: sampler
    };
}


// call this function to re-render
function drawScene(renderer) {
    let gl = renderer.gl;
    let programs = renderer.programs;
    let uniforms = renderer.uniforms;
    let buffers = renderer.buffers;
    let textures = renderer.textures;
    let framebuffers = renderer.framebuffers;

    /* RENDER TO FRAMEBUFFER */
    gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffers.ptFramebuffer.framebuffer);
    gl.viewport(0, 0, uniforms.iResolution.x, uniforms.iResolution.y);
    if (uniforms.iFrame == 0) {
        gl.clearColor(0, 1, 0.5, 0);
        gl.clearDepth(-1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }
    gl.useProgram(programs.renderProgram);

    // set position buffer
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, buffers.positionBuffer);
        var vertexPosition = gl.getAttribLocation(programs.renderProgram, 'vertexPosition');
        gl.vertexAttribPointer(vertexPosition, numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(vertexPosition);
    }

    // set shader uniforms
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iRz"), uniforms.iRz);
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iRx"), uniforms.iRx);
    gl.uniform1f(gl.getUniformLocation(programs.renderProgram, "iSc"), uniforms.iSc);
    gl.uniform2f(
        gl.getUniformLocation(programs.renderProgram, "iResolution"),
        uniforms.iResolution.x, uniforms.iResolution.y);
    gl.uniform4f(
        gl.getUniformLocation(programs.renderProgram, "iMouse"),
        uniforms.iMouse.x, uniforms.iMouse.y, uniforms.iMouse.z, uniforms.iMouse.w);
    gl.uniform1i(
        gl.getUniformLocation(programs.renderProgram, "iFrame"),
        uniforms.iFrame++);

    // set textures
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, framebuffers.ptFramebuffer.sampler);
    gl.uniform1i(gl.getUniformLocation(programs.renderProgram, "iChannel0"), 0);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, textures.treeSamplerGlass);
    gl.uniform1i(gl.getUniformLocation(programs.renderProgram, "treeSamplerGlass"), 1);
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, textures.treeSamplerContent);
    gl.uniform1i(gl.getUniformLocation(programs.renderProgram, "treeSamplerContent"), 2);

    // render to framebuffer
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }

    // copy rendered image
    gl.bindTexture(gl.TEXTURE_2D, framebuffers.ptFramebuffer.sampler);
    gl.copyTexImage2D(gl.TEXTURE_2D,
        0, gl.RGBA32F, 0, 0, uniforms.iResolution.x, uniforms.iResolution.y, 0);

    /* RENDERING TO DISPLAY */
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, uniforms.iResolution.x, uniforms.iResolution.y);
    gl.clearColor(0, 0.5, 1, 1);
    gl.clearDepth(-1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.useProgram(programs.displayProgram);

    // set position buffer
    {
        const numComponents = 2; // pull out 2 values per iteration
        const type = gl.FLOAT; // the data in the buffer is 32bit floats
        const normalize = false; // don't normalize
        const stride = 0; // how many bytes to get from one set of values to the next
        const offset = 0; // how many bytes inside the buffer to start from
        gl.bindBuffer(gl.ARRAY_BUFFER, buffers.positionBuffer);
        var vertexPosition = gl.getAttribLocation(programs.displayProgram, 'vertexPosition');
        gl.vertexAttribPointer(vertexPosition, numComponents, type, normalize, stride, offset);
        gl.enableVertexAttribArray(vertexPosition);
    }

    // set texture
    //gl.bindTexture(gl.TEXTURE_2D, framebuffers.ptFramebuffer.sampler);
    gl.bindTexture(gl.TEXTURE_2D, framebuffers.ptFramebuffer.texture);
    gl.uniform1i(gl.getUniformLocation(programs.displayProgram, "iChannel0"), 0);

    // render to canvas
    {
        const offset = 0;
        const vertexCount = 4;
        gl.drawArrays(gl.TRIANGLE_STRIP, offset, vertexCount);
    }
}


function main() {
    const canvas = document.getElementById("canvas");
    const gl = canvas.getContext("webgl2") || canvas.getContext("experimental-webgl2");
    if (gl == null) throw ("Failed to load WebGL2");

    var renderer = {
        canvas: canvas,
        gl: gl,
        extentions: {},
        programs: {
            renderProgram: null,
            displayProgram: null
        },
        uniforms: {
            iRz: -7.6 + 1e-4,
            iRx: 0.34 + 1e-4,
            iSc: 1.5,
            iResolution: { x: 1, y: 1 },
            iMouse: { x: 0, y: 0, z: 0, w: 0 },
            iFrame: 0
        },
        buffers: {
            positionBuffer: null,
        },
        textures: {
            treeSamplerGlass: null,
            treeSamplerContent: null
        },
        framebuffers: {
            ptFramebuffer: {
                framebuffer: null,
                texture: null,
                sampler: null
            }
        },
        statistics: {
            recordedTime: []
        },
        renderNeeded: true
    };
    renderer.extentions["EXT_disjoint_timer_query_webgl2"] = gl.getExtension("EXT_disjoint_timer_query_webgl2");
    if (!(renderer.extentions["EXT_color_buffer_float"] = gl.getExtension("EXT_color_buffer_float")))
        throw ("Error: your device does not support 32-bit float framebuffer.");
    renderer.uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
    renderer.uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;

    // load models first because it takes time
    loadTreeTexture(renderer, "glass1.bin", "treeSamplerGlass");
    loadTreeTexture(renderer, "content1.bin", "treeSamplerContent");

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
    //var fsSourceRender = loadShaderSource("frag-raymarch.glsl");
    var fsSourceRender = loadShaderSource("frag-octree.glsl");
    var fsSourceDisplay = `#version 300 es
        precision highp float;

        out vec4 fragColor;

        uniform sampler2D iChannel0;

        void main(void) {
            vec4 rgbn = texelFetch(iChannel0, ivec2(gl_FragCoord.xy), 0);
            fragColor = rgbn / rgbn.w;
        }
    `;
    console.timeEnd("request glsl code");

    // compile shaders
    console.time("compile shaders");
    renderer.programs.renderProgram = initShaderProgram(gl, vsSource, fsSourceRender);
    renderer.programs.displayProgram = initShaderProgram(gl, vsSource, fsSourceDisplay);
    console.timeEnd("compile shaders");

    // position buffer
    var positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    var positions = [-1, 1, 1, 1, -1, -1, 1, -1];
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(positions), gl.STATIC_DRAW);
    renderer.buffers.positionBuffer = positionBuffer;

    // framebuffer
    renderer.framebuffers.ptFramebuffer = createRenderTarget(gl, renderer.uniforms.iResolution.x, renderer.uniforms.iResolution.y);

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
                renderer.uniforms.iResolution.x + "x" + renderer.uniforms.iResolution.y + " #" + renderer.uniforms.iFrame + "<br/>" +
                (1000.0 / avrfps).toFixed(0) + " ms, " + avrfps.toFixed(1) + " fps";
            drawScene(renderer);
            //renderer.renderNeeded = false;
        }
        then = now;
        if (renderer.textures.treeSamplerGlass == null ||
            renderer.textures.treeSamplerContent == null) {
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
        renderer.uniforms.iSc *= sc;
        renderer.renderNeeded = true;
        renderer.uniforms.iFrame = 0;
    }, { passive: false });
    canvas.addEventListener("pointerdown", function (event) {
        //event.preventDefault();
        renderer.uniforms.iMouse.z = event.clientX;
        renderer.uniforms.iMouse.w = event.clientY;
        mouseDown = true;
    });
    window.addEventListener("pointerup", function (event) {
        event.preventDefault();
        renderer.uniforms.iMouse.z = -Math.abs(renderer.uniforms.iMouse.z);
        renderer.uniforms.iMouse.w = -Math.abs(renderer.uniforms.iMouse.w);
        mouseDown = false;
    });
    window.addEventListener("resize", function (event) {
        // update resolution
        renderer.uniforms.iResolution.x = canvas.width = canvas.style.width = window.innerWidth;
        renderer.uniforms.iResolution.y = canvas.height = canvas.style.height = window.innerHeight;
        // update framebuffer
        gl.deleteFramebuffer(renderer.framebuffers.ptFramebuffer.framebuffer);
        gl.deleteTexture(renderer.framebuffers.ptFramebuffer.texture);
        gl.deleteTexture(renderer.framebuffers.ptFramebuffer.sampler);
        renderer.framebuffers.ptFramebuffer = createRenderTarget(gl, renderer.uniforms.iResolution.x, renderer.uniforms.iResolution.y);
        renderer.uniforms.iFrame = 0;
    });
    canvas.addEventListener("pointermove", function (e) {
        if (mouseDown) {
            renderer.uniforms.iRx += 0.01 * e.movementY;
            renderer.uniforms.iRz -= 0.01 * e.movementX;
            renderer.uniforms.iFrame = 0;
            renderer.renderNeeded = true;
        }
        if (renderer.uniforms.iMouse.z > 0.0) {
            renderer.uniforms.iMouse.w = -Math.abs(renderer.uniforms.iMouse.w);
            renderer.uniforms.iMouse.x = e.clientX;
            renderer.uniforms.iMouse.y = e.clientY;
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