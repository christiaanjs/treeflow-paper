import pathlib
import re

import pytest

import treeflow_pipeline.submission as submission


@pytest.fixture
def manuscript_tex_file():
    path = pathlib.Path("manuscript") / "out" / "treeflow.tex"
    if not path.exists():
        pytest.skip(f"{path} has not been built")
    return path


@pytest.fixture
def manuscript_tex(manuscript_tex_file):
    return manuscript_tex_file.read_text()


@pytest.fixture
def main_text(manuscript_tex_file):
    return submission.build_main_text(manuscript_tex_file, bibliography="treeflow")


def uncommented_lines(content):
    return [line for line in content.split("\n") if not line.lstrip().startswith("%")]


def test_find_environments_skips_commented():
    content = "\n".join(
        [
            "\\begin{figure}",
            "\\caption{Real}",
            "\\end{figure}",
            "%\\begin{figure}",
            "%\\caption{Template example}",
            "%\\end{figure}",
        ]
    )
    (start, end), = submission.find_environments(content, "figure")
    assert content[start:end].endswith("\\caption{Real}\n\\end{figure}")


def test_find_environments_unterminated():
    with pytest.raises(ValueError, match="Unterminated"):
        submission.find_environments("\\begin{figure}\n", "figure")


def test_figure_legend_drops_graphics_and_widens_subfigures():
    figure = "\n".join(
        [
            "\\begin{figure}",
            "    \\centering",
            "    \\begin{subfigure}[t]{0.4\\linewidth}",
            "        \\includegraphics[width=\\linewidth]{figures/a}",
            "        \\caption{Panel a}",
            "    \\end{subfigure}",
            "    \\caption{Whole figure}",
            "    \\label{fig:example}",
            "\\end{figure}",
        ]
    )
    legend = submission.figure_legend(figure)
    assert legend.startswith("\\begin{figure}[!ht]")
    assert "includegraphics" not in legend
    assert "\\centering" not in legend
    assert "\\begin{subfigure}[t]{\\linewidth}" in legend
    assert "\\caption{Panel a}" in legend
    assert "\\caption{Whole figure}" in legend
    assert "\\label{fig:example}" in legend


def test_strip_minted_removes_package_setup():
    content = "\n".join(
        [
            "\\usepackage[cachedir=manuscript/out/minted-cache]{minted}",
            "\\usepackage{multirow}",
            "\\setminted{baselinestretch=0.9}",
            "",
        ]
    )
    assert submission.strip_minted(content) == "\\usepackage{multirow}\n"


def test_strip_minted_rejects_listings():
    content = "\\usepackage{minted}\n\\begin{minted}{python}\nx = 1\n\\end{minted}\n"
    with pytest.raises(ValueError, match="shell-escape"):
        submission.strip_minted(content)


def test_move_figures_to_legends_requires_bibliography():
    with pytest.raises(ValueError, match="bibliography"):
        submission.move_figures_to_legends(
            "\\begin{figure}\n\\caption{A}\n\\end{figure}\n"
        )


def test_main_text_has_no_figures_in_body(main_text):
    body, legends = main_text.split("\\section*{Figure Legends}")
    assert not submission.find_environments(body, "figure")
    assert "includegraphics" not in "\n".join(uncommented_lines(body))
    assert "includestandalone" not in body


def test_main_text_legends_follow_reference_list(main_text):
    assert main_text.index("\\bibliography{treeflow}") < main_text.index(
        "\\section*{Figure Legends}"
    )


def test_main_text_keeps_every_figure_legend(manuscript_tex, main_text):
    def captions(content):
        # Moving the figures past the table changes the order captions appear
        # in, so compare them as a set.
        return sorted(re.findall(r"\\caption\{(.{0,60})", content))

    def labels(content):
        return re.findall(r"\\label\{(fig:[^{}]*)\}", content)

    assert captions(main_text) == captions(manuscript_tex)
    # The figures keep their relative order, so the numbering is unchanged.
    assert labels(main_text) == labels(manuscript_tex)
    assert len(submission.find_environments(main_text, "figure")) == len(
        submission.find_environments(manuscript_tex, "figure")
    )


def test_main_text_keeps_table_embedded(manuscript_tex, main_text):
    assert len(submission.find_environments(main_text, "table")) == len(
        submission.find_environments(manuscript_tex, "table")
    )
    assert "\\begin{tabular}" in main_text


def test_main_text_needs_no_shell_escape(main_text):
    assert "minted" not in "\n".join(uncommented_lines(main_text))


def test_manifest_lists_file_names():
    content = submission.manifest(
        dict(Manuscript="out/final-submission/treeflow-main-text.tex"),
        {"figure-1": "out/final-submission/figure-1.pdf"},
        dict(Bibliography="out/final-submission/treeflow.bib"),
    )
    assert "treeflow-main-text.tex: Manuscript" in content
    assert "figure-1.pdf: figure-1" in content
    assert "treeflow.bib: Bibliography" in content
