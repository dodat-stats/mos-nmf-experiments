#!/usr/bin/env python3

"""Combine the phone-readable MiSo PNG pages into one PDF."""

from pathlib import Path

from reportlab.lib.pagesizes import portrait
from reportlab.pdfgen import canvas
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "output" / "sci-space-miso-large" / "phone-report"
PAGES = [
    OUTPUT / "01-goals-summary.png",
    OUTPUT / "02-factor-graph.png",
    OUTPUT / "03-factor-interpretation-a.png",
    OUTPUT / "04-factor-interpretation-b.png",
]
PDF = OUTPUT / "miso-sci-space-findings-phone.pdf"


def main() -> None:
    # A 6 x 8.67 inch portrait page matches the PNG aspect ratio and reads
    # comfortably in a phone PDF viewer.
    page_size = portrait((432.0, 624.0))
    document = canvas.Canvas(str(PDF), pagesize=page_size)
    document.setTitle("MiSo sci-Space findings and factor graph")
    document.setAuthor("MiSo project")
    document.setSubject("Goal evaluation and interpretation of MiSo factors")

    page_width, page_height = page_size
    for page_path in PAGES:
        if not page_path.exists():
            raise FileNotFoundError(page_path)
        with Image.open(page_path) as image:
            width, height = image.size
        scale = min(page_width / width, page_height / height)
        draw_width = width * scale
        draw_height = height * scale
        x = (page_width - draw_width) / 2
        y = (page_height - draw_height) / 2
        document.drawImage(
            str(page_path),
            x,
            y,
            width=draw_width,
            height=draw_height,
            preserveAspectRatio=True,
            mask="auto",
        )
        document.showPage()

    document.save()
    print(PDF)


if __name__ == "__main__":
    main()

