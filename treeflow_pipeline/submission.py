"""Helpers for assembling the final journal submission from the built manuscript.

Systematic Biology asks for the main text as a clean, editable LaTeX source
(accompanied by its PDF, ``.bib`` and class files), with the figures removed
from the text and uploaded as separate files, and the figure legends collected
after the reference list.  Tables may stay embedded as long as they remain
editable, which they are: ``build_manuscript`` inlines the generated table
source into the manuscript.

The transformations below operate on the compiled ``manuscript/out/treeflow.tex``
rather than re-rendering the Jinja templates, so the submitted source is the
same source that produced the manuscript PDF.  The figure ``figure``
environments are moved rather than rewritten into a plain list, which keeps the
``\\label``/``\\ref`` machinery (including subfigure letters) working and keeps
the figure numbering identical to the compiled manuscript.
"""
import pathlib
import re

# Placement for the moved figure floats.  ``!ht`` asks LaTeX to typeset each
# legend where it appears in the legends section rather than deferring it to a
# float page, which keeps them in figure order under the section heading.
LEGEND_PLACEMENT = "[!ht]"
LEGENDS_HEADING = "Figure Legends"

GRAPHICS_COMMANDS = ("includegraphics", "includestandalone")

_GRAPHICS_LINE = re.compile(
    r"^[ \t]*\\(?:" + "|".join(GRAPHICS_COMMANDS) + r")(?:\[[^\]]*\])?\{[^{}]*\}[ \t]*\n?",
    re.MULTILINE,
)
_CENTERING_LINE = re.compile(r"^[ \t]*\\centering[ \t]*\n", re.MULTILINE)
_SUBFIGURE_WIDTH = re.compile(r"(\\begin\{subfigure\}(?:\[[^\]]*\])?)\{[^{}]*\}")
_FIGURE_PLACEMENT = re.compile(r"\A\\begin\{figure\}(?:\[[^\]]*\])?")
_BIBLIOGRAPHY = re.compile(r"^\\bibliography\{[^{}]*\}[ \t]*$", re.MULTILINE)
_MINTED_PREAMBLE = re.compile(
    r"^\\(?:usepackage(?:\[[^\]]*\])?\{minted\}"
    r"|setminted(?:\[[^\]]*\])?\{[^{}]*\})[ \t]*\n",
    re.MULTILINE,
)


def _is_commented(content, position):
    """Is ``position`` preceded on its line by an unescaped comment character?"""
    line_start = content.rfind("\n", 0, position) + 1
    return re.search(r"(?<!\\)%", content[line_start:position]) is not None


def find_environments(content, name):
    """Spans of every uncommented ``name`` environment, as ``(start, end)`` pairs.

    Environments of the same name are assumed not to nest, which holds for the
    manuscript's ``figure`` (the nested environments are ``subfigure``) and
    ``minted`` environments.
    """
    begin = re.compile(r"\\begin\{" + name + r"\}")
    end = re.compile(r"\\end\{" + name + r"\}")
    spans = []
    position = 0
    while True:
        begin_match = begin.search(content, position)
        if begin_match is None:
            return spans
        if _is_commented(content, begin_match.start()):
            position = begin_match.end()
            continue
        end_match = end.search(content, begin_match.end())
        if end_match is None:
            raise ValueError(f"Unterminated {name} environment in manuscript source")
        spans.append((begin_match.start(), end_match.end()))
        position = end_match.end()


def _strip_graphics(figure):
    """Remove the image inclusion commands from a figure environment."""
    return _GRAPHICS_LINE.sub("", figure)


def _widen_subfigures(figure):
    """Give every subfigure the full text width.

    Subfigure boxes sized for images side by side (``0.4\\linewidth``) leave the
    subcaptions in narrow columns once the images are gone; at full width they
    read as consecutive paragraphs of the legend.
    """
    return _SUBFIGURE_WIDTH.sub(r"\1{\\linewidth}", figure)


def _set_placement(figure, placement):
    return _FIGURE_PLACEMENT.sub(r"\\begin{figure}" + placement, figure, count=1)


def figure_legend(figure, placement=LEGEND_PLACEMENT):
    """Turn a figure environment into its caption-only legend.

    The ``\\centering`` that positioned the (now absent) image is dropped along
    with it, leaving a float that holds nothing but the caption.
    """
    legend = _CENTERING_LINE.sub("", _strip_graphics(figure))
    return _set_placement(_widen_subfigures(legend), placement)


def move_figures_to_legends(
    content, heading=LEGENDS_HEADING, placement=LEGEND_PLACEMENT
):
    """Remove the figures from the text and list their legends after the references.

    Each ``figure`` environment is stripped of its images and moved, in order,
    into a figure legends section following the ``\\bibliography`` command.

    Parameters
    ----------
    content : str
        Contents of the compiled manuscript ``.tex`` file.
    heading : str
        Title of the section holding the legends.
    placement : str
        Float placement specifier applied to the moved figures.

    Returns
    -------
    str
        The manuscript source with figures replaced by a legends section.
    """
    spans = find_environments(content, "figure")
    if not spans:
        raise ValueError("No figure environments found in manuscript source")
    legends = [figure_legend(content[start:end], placement) for start, end in spans]

    body = content
    for start, end in reversed(spans):
        # Swallow the newline that terminated the figure environment so the
        # blank lines around it collapse into a single paragraph break.
        stop = end + 1 if body[end : end + 1] == "\n" else end
        body = body[:start] + body[stop:]

    bibliography = _BIBLIOGRAPHY.search(body)
    if bibliography is None:
        raise ValueError("No \\bibliography command found in manuscript source")
    legends_section = "\n".join(
        ["", "", "\\clearpage", "", f"\\section*{{{heading}}}", ""]
        + [legend + "\n" for legend in legends]
    )
    return body[: bibliography.end()] + legends_section + body[bibliography.end() :]


def strip_minted(content):
    """Drop the ``minted`` setup from the submission copy of the manuscript.

    ``minted`` requires ``pdflatex --shell-escape``, which a publisher's
    production system cannot be expected to run.  The main text no longer
    contains code listings, so the package is simply removed; if a listing is
    reintroduced this raises rather than silently emitting a source that the
    journal cannot compile.
    """
    minted_environments = find_environments(content, "minted")
    if minted_environments:
        raise ValueError(
            "Manuscript contains minted environments, which need "
            "pdflatex --shell-escape; freeze the minted cache for the "
            "submission copy instead of stripping the package"
        )
    return _MINTED_PREAMBLE.sub("", content)


def set_bibliography(content, bibliography):
    """Point ``\\bibliography`` at ``bibliography`` (a ``.bib`` file stem)."""
    if not _BIBLIOGRAPHY.search(content):
        raise ValueError("No \\bibliography command found in manuscript source")
    return _BIBLIOGRAPHY.sub(f"\\\\bibliography{{{bibliography}}}", content)


def build_main_text(tex_file, bibliography=None):
    """Build the clean main text source for submission.

    Removes the figures (uploaded separately) and lists their legends after the
    reference list, and removes the shell-escape-dependent ``minted`` setup.

    Parameters
    ----------
    tex_file : str or pathlib.Path
        The compiled manuscript ``.tex`` file.
    bibliography : str, optional
        Stem of the ``.bib`` file as it is named in the submission bundle.

    Returns
    -------
    str
        Contents of the submission main text ``.tex`` file.
    """
    content = pathlib.Path(str(tex_file)).read_text()
    content = strip_minted(content)
    content = move_figures_to_legends(content)
    if bibliography is not None:
        content = set_bibliography(content, bibliography)
    return content


def manifest(main_text, figures, supporting):
    """Describe the submission bundle, one line per file.

    Parameters
    ----------
    main_text : dict
        Mapping of description to path for the main text files.
    figures : dict
        Mapping of description to path for the separate figure files.
    supporting : dict
        Mapping of description to path for the remaining files.

    Returns
    -------
    str
        Contents of the bundle's manifest file.
    """

    def section(title, entries):
        return [title, "-" * len(title)] + [
            f"{pathlib.Path(str(path)).name}: {description}"
            for description, path in entries.items()
        ]

    return (
        "\n".join(
            ["Final submission files", "======================", ""]
            + section("Main text", main_text)
            + [""]
            + section("Figures (one file per figure)", figures)
            + [""]
            + section("Supporting files", supporting)
        )
        + "\n"
    )
