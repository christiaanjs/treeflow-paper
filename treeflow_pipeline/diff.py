"""Helper functions for the ms.smk Snakemake workflow."""
import re
import subprocess


LATEXDIFF_ARGS = [
    "--append-safecmd=bibliographystyle,setcitestyle,usetikzlibrary,definecolor",
    "--append-context2cmd=maketitle",
    "--add-to-config", "VERBATIMENV=minted",
]


def _fix_deleted_verbatim(content):
    """Comment out verbatim code lines that leak out of deleted minted blocks.

    When a figure containing a minted environment is entirely deleted,
    latexdiff wraps the figure in DIFdelbegin/DIFdelend but the verbatim
    content of the minted block is not automatically commented out,
    causing LaTeX errors (e.g. bare underscores outside math mode).
    """
    lines = content.split("\n")
    out, in_del, in_verb = [], False, False
    for line in lines:
        s = line.rstrip()
        if r"\DIFdelbegin" in s and not s.startswith("%"):
            in_del = True
        if r"\DIFdelend" in s and not s.startswith("%"):
            in_del = False
            in_verb = False
        if in_del and re.match(r"%DIFDELCMD\s*<\s*\\begin\{minted", s):
            in_verb = True
        if in_del and re.match(r"%DIFDELCMD\s*<\s*\\end\{minted", s):
            in_verb = False
        if in_del and in_verb and not s.startswith("%") and s:
            out.append("%DIFDELCMD < " + line)
        else:
            out.append(line)
    return "\n".join(out)


def run_latexdiff(old_tex, new_tex):
    """Run latexdiff and return post-processed diff content.

    Handles three incompatibilities between the old article-class submission
    and the new sysbio_sse journal-class submission:

    1. Removes the ``DELETED TITLE COMMANDS FOR MARKUP`` preamble block,
       which conflicts with sysbio_sse's ``\\author{}`` using ``\\@dblarg``.
    2. Replaces the body title/header section with the clean new version,
       because ``\\DIFadd{}`` wrappers inside ``\\author{}`` and around
       ``\\maketitle`` break the journal class.
    3. Comments out verbatim code inside deleted figures and removes
       entirely-deleted figure environments that have brace-balance issues
       in their captions.

    Parameters
    ----------
    old_tex : str
        Path to the old .tex file (base for comparison).
    new_tex : str
        Path to the new .tex file (target version).

    Returns
    -------
    str
        Post-processed diff .tex content ready for pdflatex compilation.
    """
    result = subprocess.run(
        ["latexdiff"] + LATEXDIFF_ARGS + [str(old_tex), str(new_tex)],
        capture_output=True, text=True, check=True,
    )
    raw = result.stdout

    raw = _fix_deleted_verbatim(raw)

    # Remove entirely-deleted figure environments (brace-balance issues in captions)
    raw = re.sub(
        r"\\DIFdelbegin %DIFDELCMD < \\begin\{figure\}.*?\\DIFdelend\s*\n",
        "",
        raw,
        flags=re.DOTALL,
    )

    # Remove "DELETED TITLE COMMANDS FOR MARKUP" block from preamble
    raw = re.sub(
        r"%DIF DELETED TITLE COMMANDS FOR MARKUP\n.*?(?=%DIF PREAMBLE EXTENSION)",
        "",
        raw,
        flags=re.DOTALL,
    )

    # Replace the body title/header with the clean new version
    with open(str(new_tex)) as fh:
        new_tex_content = fh.read()
    new_header_m = re.search(
        r"(\\begin\{document\}.*?(?=\\section\{Introduction\}))",
        new_tex_content, re.DOTALL,
    )
    diff_body_m = re.search(r"(\\section\{Introduction\}.*)", raw, re.DOTALL)
    preamble_m = re.search(r"(.*?)(?=\\begin\{document\})", raw, re.DOTALL)
    return preamble_m.group(1) + new_header_m.group(1) + diff_body_m.group(1)


def extract_table_tex(treeflow_tex, remove_caption=False):
    """Extract the first table environment from treeflow.tex as a standalone document.

    Parameters
    ----------
    treeflow_tex : str
        Path to the compiled treeflow.tex manuscript file.
    remove_caption : bool
        If True, strip ``\\caption{...}`` and ``\\label{...}`` from the
        extracted table body before wrapping.

    Returns
    -------
    str
        Complete standalone LaTeX document wrapping the extracted table.
    """
    with open(str(treeflow_tex)) as fh:
        content = fh.read()
    m = re.search(r"(\\begin\{table\}.*?\\end\{table\})", content, re.DOTALL)
    if not m:
        raise ValueError(f"No table environment found in {treeflow_tex}")
    table_body = m.group(1)
    if remove_caption:
        table_body = re.sub(
            r"\s*\\caption\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}", "", table_body
        )
        table_body = re.sub(r"\s*\\label\{[^{}]*\}", "", table_body)
    return (
        "\\documentclass{article}\n"
        "\\usepackage{booktabs}\n"
        "\\usepackage{multirow}\n"
        "\\begin{document}\n"
        + table_body + "\n"
        "\\end{document}\n"
    )


def extract_compound_figure_tex(treeflow_tex, figure_index):
    """Extract a compound figure (containing subfigures) as a standalone document.

    Includes the subfigure images and their subcaptions but omits the main
    figure-level caption, as required for separate figure file submission.

    Parameters
    ----------
    treeflow_tex : str
        Path to the compiled treeflow.tex manuscript file.
    figure_index : int
        Zero-based index among figures that contain subfigures.

    Returns
    -------
    str
        Complete standalone LaTeX document wrapping the extracted figure.
    """
    with open(str(treeflow_tex)) as fh:
        content = fh.read()

    figure_envs = re.findall(
        r"\\begin\{figure\}.*?\\end\{figure\}", content, re.DOTALL
    )
    compound = [f for f in figure_envs if r"\begin{subfigure}" in f]
    if figure_index >= len(compound):
        raise ValueError(
            f"figure_index {figure_index} out of range "
            f"({len(compound)} compound figures found)"
        )
    body = compound[figure_index]

    # Remove the main \caption{...} that follows the last \end{subfigure}.
    # Handles one level of nested braces (sufficient for plain-text captions).
    last_sub_end = body.rfind(r"\end{subfigure}") + len(r"\end{subfigure}")
    tail = re.sub(
        r"\\caption\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}",
        "",
        body[last_sub_end:],
    )
    body = body[:last_sub_end] + tail

    return (
        "\\documentclass{article}\n"
        "\\usepackage{graphicx}\n"
        "\\usepackage{subcaption}\n"
        "\\begin{document}\n"
        + body + "\n"
        "\\end{document}\n"
    )

