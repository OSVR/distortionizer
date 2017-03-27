#!/bin/sh

# Takes each line as being a filename, turns it into a null-delimited string
# and hands it to xargs -0 to pass to clang-format efficiently (because running
# clang-format once for each file is very slow)
run_clang_format_on_input_lines() {
    (while read fn; do
        printf "%s\0" "$fn"
    done) | xargs -0 clang-format -style=file -i
}

(
cd $(dirname $0)

for dir in angles_to_config calibration shaders; do
(
    cd $dir
    ls *.h *.cpp | run_clang_format_on_input_lines
)
done

for dir in RenderManager_image_based RenderManager_xbox; do
(
    cd RenderManager_image_based
    ls *.cpp | run_clang_format_on_input_lines
)
done

)
