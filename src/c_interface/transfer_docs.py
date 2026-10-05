#!/usr/bin/env python3
"""Transfer documentation from Fortran source files to C interface files."""

import os
import re
import glob

SRC_DIR = "src"
C_INTERFACE_DIR = os.path.join(SRC_DIR, "c_interface")


def extract_fortran_doc(fortran_path):
    """Extract the doxygen documentation block from a Fortran file.

    Returns the documentation lines (including !> and ! prefixes) from
    the first line matching '!> \\brief' up to (but not including) the
    'subroutine' statement line.
    """
    with open(fortran_path, 'r') as f:
        lines = f.readlines()

    doc_lines = []
    in_doc = False
    for line in lines:
        stripped = line.rstrip('\n')
        if not in_doc:
            # Look for the start of doxygen documentation
            if re.match(r'^!>\s*\\brief', stripped):
                in_doc = True
                doc_lines.append(stripped)
        else:
            # Stop at the subroutine statement
            if re.match(r'^subroutine\s', stripped):
                break
            doc_lines.append(stripped)

    return doc_lines


def lowercase_param_names(content):
    """Convert parameter names in \\param lines to lowercase.

    \\param[in] R  -> \\param[in] r
    \\param[out] Q -> \\param[out] q
    """
    return re.sub(
        r'(\\param\[[^\]]*\]\s+)(\w+)',
        lambda m: m.group(1) + m.group(2).lower(),
        content
    )


def convert_group_name(content):
    """Convert Fortran group names to C group names for \\ingroup tags.

    \\ingroup qrdecomp    -> \\ingroup c_qrdecomp
    \\ingroup choldecomp  -> \\ingroup c_choldecomp
    \\ingroup ludecomp    -> \\ingroup c_ludecomp
    """
    group_mapping = {
        'qrdecomp': 'c_qrdecomp',
        'choldecomp': 'c_choldecomp',
        'ludecomp': 'c_ludecomp',
    }

    def replace_group(m):
        group_name = m.group(1)
        new_name = group_mapping.get(group_name, group_name)
        return f'\\ingroup {new_name}'

    return re.sub(r'\\ingroup\s+(\w+)', replace_group, content)


def extract_c_signature(c_path):
    """Extract the C function signature from a C file.

    Returns the signature string (e.g., 'void qrupdate_caxcpy(...)') or None if not found.
    """
    with open(c_path, 'r') as f:
        content = f.read()

    # Find the QRUPDATE_EXPORT function signature
    # Pattern: QRUPDATE_EXPORT void funcname(...)
    # The signature may span multiple lines
    match = re.search(r'QRUPDATE_EXPORT\s+void\s+\w+\s*\([^)]*\)', content, re.DOTALL)
    if match:
        signature = match.group(0)
        # Clean up whitespace: collapse multiple spaces/newlines into single spaces
        signature = re.sub(r'\s+', ' ', signature).strip()
        return signature
    return None


def convert_to_c_doc(fortran_doc_lines, c_signature=None):
    """Convert Fortran doxygen comment lines to C doxygen comment block.

    !> text  -> * text
    ! text   -> * text   (section underlines like ! =============)
    !>       -> *        (empty doxygen lines)
    !        -> *        (empty Fortran comment lines within doc)

    Parameter names in \\param lines are converted to lowercase.
    If c_signature is provided, the \\par Definition block is replaced
    with the C function signature.
    """
    c_lines = ["/**"]
    in_definition_block = False
    definition_block_ended = False
    for line in fortran_doc_lines:
        # Remove the trailing newline if present
        # Match patterns:
        # "!> text" -> "* text"
        # "!> " (empty) -> "*"
        # "! text" -> "* text"
        # "!" (just exclamation) -> "*"
        m = re.match(r'^!>\s?(.*)', line)
        if m:
            content = m.group(1)
        else:
            m = re.match(r'^!\s?(.*)', line)
            if m:
                content = m.group(1)
            else:
                continue

        # Check for start of definition block
        if c_signature and re.match(r'\\par Definition:', content):
            in_definition_block = True
            # Insert the C signature in a verbatim block
            c_lines.append("  \\par C Interface:")
            c_lines.append("  ==============")
            c_lines.append("  \\verbatim")
            c_lines.append(f"    {c_signature}")
            c_lines.append("  \\endverbatim")
            continue

        # Skip lines inside the definition block
        if in_definition_block:
            if re.match(r'\\endverbatim', content):
                in_definition_block = False
                definition_block_ended = True
            continue

        if content:
            content = lowercase_param_names(content)
            content = convert_group_name(content)
            c_lines.append(f"  {content}")
        else:
            c_lines.append("")

    c_lines.append(" */")
    return c_lines


def strip_existing_doc(lines):
    """Strip all existing /** ... */ documentation blocks before QRUPDATE_EXPORT.

    Returns the lines with all doc blocks removed.
    """
    # Repeatedly strip doc blocks until none remain
    while True:
        # Find the QRUPDATE_EXPORT line (the real one, not indented/verbatim content)
        export_idx = None
        for i, line in enumerate(lines):
            # Match actual QRUPDATE_EXPORT declarations at column 0
            if re.match(r'^QRUPDATE_EXPORT\s+void\b', line):
                export_idx = i
                break

        if export_idx is None:
            break

        # Find the last */ before export_idx - that's the end of the last doc block
        last_end = None
        for i in range(export_idx - 1, -1, -1):
            if lines[i].strip() == '*/':
                last_end = i
                break

        if last_end is None:
            break

        # Find the corresponding /** - scan backwards looking for /** at column 0
        # but skip any /** inside \verbatim blocks
        start_doc = None
        in_verbatim = False
        for i in range(last_end - 1, -1, -1):
            stripped = lines[i].strip()
            if stripped == '\\endverbatim':
                in_verbatim = True
            elif stripped == '\\verbatim':
                in_verbatim = False
            elif stripped == '/**' and not in_verbatim:
                start_doc = i
                break

        if start_doc is None:
            break

        # Remove lines from start_doc to last_end (inclusive)
        # Ensure there's exactly one blank line before the doc block
        if start_doc > 0 and lines[start_doc - 1].strip() == '':
            lines = lines[:start_doc] + lines[export_idx:]
        else:
            lines = lines[:start_doc] + ['\n'] + lines[export_idx:]

    return lines


def insert_doc_in_c_file(c_path, c_doc_lines):
    """Insert the C documentation block before the QRUPDATE_EXPORT function."""
    with open(c_path, 'r') as f:
        lines = f.readlines()

    # Strip any existing documentation block
    lines = strip_existing_doc(lines)

    # Find the line with QRUPDATE_EXPORT
    insert_idx = None
    for i, line in enumerate(lines):
        if 'QRUPDATE_EXPORT' in line:
            insert_idx = i
            break

    if insert_idx is None:
        print(f"  WARNING: No QRUPDATE_EXPORT found in {c_path}")
        return False

    # Build the new content
    new_lines = lines[:insert_idx]
    for doc_line in c_doc_lines:
        new_lines.append(doc_line + '\n')
    new_lines.append('\n')
    new_lines.extend(lines[insert_idx:])

    with open(c_path, 'w') as f:
        f.writelines(new_lines)

    return True


def main():
    # Find all Fortran files matching [cdsz]*.f90
    fortran_pattern = os.path.join(SRC_DIR, "[cdsz]*.f90")
    fortran_files = sorted(glob.glob(fortran_pattern))

    processed = 0
    skipped = 0

    for f90_path in fortran_files:
        basename = os.path.splitext(os.path.basename(f90_path))[0]
        c_path = os.path.join(C_INTERFACE_DIR, f"{basename}.c")

        if not os.path.exists(c_path):
            print(f"SKIP {basename}: no C interface file found")
            skipped += 1
            continue

        # Extract Fortran documentation
        doc_lines = extract_fortran_doc(f90_path)
        if not doc_lines:
            print(f"SKIP {basename}: no documentation found in Fortran file")
            skipped += 1
            continue

        # Extract C function signature
        c_signature = extract_c_signature(c_path)

        # Convert to C documentation
        c_doc = convert_to_c_doc(doc_lines, c_signature)

        # Insert in C file
        if insert_doc_in_c_file(c_path, c_doc):
            print(f"DONE {basename}")
            processed += 1

    print(f"\nProcessed: {processed}, Skipped: {skipped}")


if __name__ == '__main__':
    main()
