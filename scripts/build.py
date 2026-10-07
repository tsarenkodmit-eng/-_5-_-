#!/usr/bin/env python3
"""Build each theory part as an independent PDF; never concatenate the sources."""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PARTS = {
    1: "theory-01-linear-algebra",
    2: "theory-02-metric-curvature",
    3: "theory-03-constant-curvature-symmetry",
}


def to_standard_markdown(text: str) -> str:
    text = re.sub(r"\$`([^`\n]+)`\$", r"$\1$", text)
    return re.sub(r"^```math\n(.*?)\n```$", r"$$\n\1\n$$", text,
                  flags=re.MULTILINE | re.DOTALL)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", type=int, choices=PARTS,
                        help="Build only this part (default: all parts separately)")
    parser.add_argument("--obsidian", type=Path,
                        help="Export one part with standard math delimiters; requires --part")
    args = parser.parse_args()
    if args.obsidian and args.part is None:
        parser.error("--obsidian requires --part; parts are exported separately")
    selected = [args.part] if args.part is not None else list(PARTS)
    for number in selected:
        stem = PARTS[number]
        text = to_standard_markdown((ROOT / f"{stem}.md").read_text(encoding="utf-8"))
        if args.obsidian:
            args.obsidian.write_text(text, encoding="utf-8")
            print(f"Saved: {args.obsidian.resolve()}")
            continue
        title, text = text.split("\n", 1)
        subtitle = title.removeprefix("# ")
        with tempfile.TemporaryDirectory(prefix=f"geometry-part-{number}-") as temp:
            source = Path(temp) / "part.md"
            source.write_text(text, encoding="utf-8")
            subprocess.run([
                "pandoc", str(source), "--from=markdown+tex_math_dollars",
                "--pdf-engine=xelatex", "--metadata-file", str(ROOT / "config/pdf.yaml"),
                "--metadata", f"subtitle={subtitle}",
                "--include-in-header", str(ROOT / "config/header.tex"),
                "--table-of-contents", "--toc-depth=2", "--standalone",
                "--output", str(ROOT / f"{stem}.pdf"),
            ], check=True, cwd=ROOT)
        print(f"Built: {ROOT / (stem + '.pdf')}")


if __name__ == "__main__":
    main()
