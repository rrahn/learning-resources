# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (lr_test_snippet CXX)

include (../ci.cmake)

macro (lr_snippet test_name_prefix snippet snippet_base_path)
    lr_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    lr_test_component (snippet_test_name "${snippet}" TEST_NAME)
    lr_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} lr::test::unit)
    set_target_properties(${target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${snippet_target_path}"
    )

    add_test (NAME "${test_name_prefix}/${snippet_test_name}_snippet" COMMAND ${target})

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (lr_snippets test_name_prefix snippet_base_path)
    lr_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        lr_snippet("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

lr_require_ccache ()
lr_require_test ()

lr_snippets ("snippet" "${CMAKE_SOURCE_DIR}")
