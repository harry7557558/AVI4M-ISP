# Convert GLSL code to C++ code so they can be compiled after `#include "glsl2cpp.h"`


import sys
import re
import os


def convert_swizzling(glsl_source: str) -> str:

    # get a list of possible swizzlings
    sws = []
    comps = ['x', 'y', 'z', 'w']  # please only use x/y/z/w and not r/g/b/s/t/p
    for l in range(4, 1, -1):
        for i in range(4**l):
            swizzle = ""
            for b in range(l):
                swizzle += comps[(i//(4**b)) % 4]
            sws.append('.'+swizzle[::-1])
    assert len(sws) == 256+64+16

    # replace swizzlings
    for sw in sws:
        if sw not in glsl_source:
            continue
        regex = re.compile(f"(\\{sw})([^0-9A-Za-z_\(\[])")
        glsl_source = regex.sub("\\1()\\2", glsl_source)

    return glsl_source


def replace_inout(glsl_source: str) -> str:

    # remove `in`
    regex = re.compile(
        "([\(\,\s])in\s+([A-Za-z0-9_]+\s+[A-Za-z0-9_]+\s*[\,\)])")
    glsl_source = regex.sub("\\1\\2", glsl_source)
    glsl_source = regex.sub("\\1\\2", glsl_source)

    # replace `out` and `inout`
    regex = re.compile(
        "([\(\,\s])(out|inout)\s+([A-Za-z0-9_]+)\s+([A-Za-z0-9_]+\s*[\,\)])")
    glsl_source = regex.sub("\\1\\3 &\\4", glsl_source)
    glsl_source = regex.sub("\\1\\3 &\\4", glsl_source)

    return glsl_source


def convert_float(glsl_source: str) -> str:

    # convert #.# to #.#f
    regex = re.compile("([^A-Za-z0-9_])([0-9]*\\.[0-9]+)([^A-Za-z0-9_])")
    glsl_source = regex.sub("\\1\\2f\\3", glsl_source)
    glsl_source = regex.sub("\\1\\2f\\3", glsl_source)
    regex = re.compile("([^A-Za-z0-9_])([0-9]+\\.[0-9]*)([^A-Za-z0-9_])")
    glsl_source = regex.sub("\\1\\2f\\3", glsl_source)
    glsl_source = regex.sub("\\1\\2f\\3", glsl_source)

    # scientific notation: ??

    return glsl_source


def glsl2cpp(filepath: str) -> str:
    directory = filepath[:filepath.rfind('/')+1]
    print("Open GLSL source:", filepath)
    glsl_source = open(filepath, "r").read()
    glsl_source = glsl_source.replace('\t', ' '*4)
    glsl_source = convert_swizzling(glsl_source)
    glsl_source = replace_inout(glsl_source)
    glsl_source = convert_float(glsl_source)

    glsl_source = glsl_source.replace('\r\n', '\n').split('\n')
    glsl_source_1 = []
    for line in glsl_source:
        if line.startswith('#include'):
            # assume relative path
            regex = r"#include [\<\"\']([A-Za-z0-9\-_\.\/\\]+)[\>\"\']"
            sourcepath = directory + re.match(regex, line).group(1)
            glsl_source_1.append(glsl2cpp(sourcepath))
        elif line.startswith('#iChannel'):
            regex = r"#(iChannel\d) [\<\"\']([A-Za-z0-9\-_\.\/\\]+)[\>\"\']"
            match = re.match(regex, line)
            sampler = match.group(1)
            sourcepath = os.path.abspath(directory + match.group(2)).replace('\\', '/')
            glsl_source_1.append(f"const sampler2D {sampler}(\"{sourcepath}\");")
        else:
            glsl_source_1.append(line)
    glsl_source = '\n'.join(glsl_source_1)

    return glsl_source.strip()


if __name__ == "__main__":
    argv = sys.argv
    if len(argv) < 2:
        argv = [__file__, ".glsl", ".glsl.cpp"]

    source = glsl2cpp(argv[1].replace('\\', '/'))

    print("Write C++ source to:", argv[2])
    with open(argv[2], "w") as fp:
        fp.write(source)
