#!/usr/bin/env python3
"""Build the checked-in PDF from the GitHub-friendly Markdown source."""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def to_standard_markdown(text: str) -> str:
    return re.sub(r"^```math\n(.*?)\n```$", r"$$\n\1\n$$", text,
                  flags=re.MULTILINE | re.DOTALL)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obsidian", type=Path,
                        help="Export an optional Markdown copy with $$ math blocks")
    args = parser.parse_args()
    text = to_standard_markdown((ROOT / "теория.md").read_text(encoding="utf-8"))
    if args.obsidian:
        args.obsidian.write_text(text, encoding="utf-8")
        print(f"Saved: {args.obsidian.resolve()}")
        return
    # The PDF title is supplied separately; all remaining heading levels move up one.
    text = re.sub(r"\A# [^\n]+\n", "", text)
    text = re.sub(r"^(#{2,}) ", lambda m: m[1][1:] + " ", text, flags=re.MULTILINE)
    text = re.sub(r"^(# Часть II\.[^\n]+)$", r"\\clearpage\n\n\1", text, flags=re.MULTILINE)
    with tempfile.TemporaryDirectory(prefix="geometry-build-") as temp:
        source = Path(temp) / "theory.md"
        source.write_text(text, encoding="utf-8")
        subprocess.run([
            "pandoc", str(source), "--from=markdown+tex_math_dollars",
            "--pdf-engine=xelatex", "--metadata-file", str(ROOT / "config/pdf.yaml"),
            "--include-in-header", str(ROOT / "config/header.tex"),
            "--table-of-contents", "--toc-depth=2", "--standalone",
            "--output", str(ROOT / "теория.pdf"),
        ], check=True, cwd=ROOT)
    print(f"Built: {ROOT / 'теория.pdf'}")


if __name__ == "__main__":
    main()
