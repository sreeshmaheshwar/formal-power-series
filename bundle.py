#!/usr/bin/env python3

"""
bundle.py

In competitive programming, a single, self-contained source
file is typically submitted to the judge. So, we provide this
makeshift bundling tool that expands THIS library's headers
included by a source file, similar to (and delegating to)
`expander.py` of Atcoder's ac-library, which can be found at
(https://github.com/atcoder/ac-library/blob/master/expander.py).

Outputs bundled, submission-ready source code to `combined.cpp`.
"""

import subprocess
import sys


FPS_H = "FormalPowerSeries.h"
FPS_CPP = "FormalPowerSeries.cpp"


def is_include_line(line, headerName):
    # FIXME: This is a very naive check. Lines may be comments, the file name
    # might be a substring of another file, etc.
    return line.startswith("#include") and headerName in line.split()[1]


def bundle_with_function_expansion(filePath, headerName, expansionFunction):
    with open(filePath) as f:
        expanded = False
        output = []

        for line in f.readlines():
            line = line.rstrip()
            if not expanded and is_include_line(line, headerName):
                expanded = True
                output.extend(expansionFunction())
            else:
                output.append(line)

        return output


def get_expanded_fps_lib(fpsLibDir):
    return bundle_with_function_expansion(
        f"{fpsLibDir}/{FPS_H}",
        FPS_CPP,
        lambda: bundle_with_function_expansion(
            f"{fpsLibDir}/{FPS_CPP}", FPS_H, lambda: []
        ),  # Remove include of FPS_H in FPS_CPP if present.
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <source_file> <fps_lib_dir>")
        sys.exit(1)

    filePath, fpsLibDir = sys.argv[1], sys.argv[2]

    with open("combined.cpp", "w") as of:
        of.write(
            "\n".join(
                bundle_with_function_expansion(
                    filePath, FPS_H, lambda: get_expanded_fps_lib(fpsLibDir)
                )
            )
        )
        subprocess.run(
            [
                "python3",
                f"{fpsLibDir}/ac-library/expander.py",
                f"--lib={fpsLibDir}/ac-library",
                "combined.cpp",
            ]
        )
