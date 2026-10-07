#!/usr/bin/env python3
"""Typeset the photographed assignment; preserve text and its original hierarchy."""
from pathlib import Path
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.colors import HexColor

ROOT = Path(__file__).resolve().parents[1]
FONT_DIR = Path('/usr/share/fonts/truetype/dejavu')
for name, filename in [('AssignmentSans','DejaVuSans.ttf'),
                       ('AssignmentSerif','DejaVuSerif.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONT_DIR / filename)))


def justified(c, text, x, y, width, font, size):
    words = text.split(' ')
    text_width = sum(pdfmetrics.stringWidth(w, font, size) for w in words)
    gap = (width-text_width)/(len(words)-1)
    c.setFont(font,size)
    for word in words:
        c.drawString(x,y,word)
        x += pdfmetrics.stringWidth(word,font,size)+gap


def main():
    destination=ROOT/'assignment-09.pdf'
    w,h=landscape(A4)
    c=canvas.Canvas(str(destination),pagesize=(w,h))
    c.setTitle('Курсовая работа — вариант №9')
    c.setAuthor('')
    c.setFont('AssignmentSans',21)
    c.drawCentredString(w/2,h-57,'КУРСОВАЯ РАБОТА')
    c.setFont('AssignmentSans',16.8)
    c.drawCentredString(w/2,h-106,'Дифференциальная геометрия и основы тензорного исчисления')
    c.setFont('AssignmentSans',21)
    c.drawCentredString(w/2,h-166,'ВАРИАНТ №9.')
    left=52;right=w-66
    c.setFont('AssignmentSans',19)
    c.drawString(left,h-215,'ТЕОРИЯ.')
    justified(c,'Псевдоевклидово пространство, определение и основные',left,h-254,
              right-left,'AssignmentSans',18)
    c.setFont('AssignmentSans',18)
    c.drawString(left,h-289,'свойства. Связь со специальной теорией относительности.')
    c.setFont('AssignmentSans',19)
    c.drawString(left,h-355,'ЗАДАНИЕ.')
    box_x=87;box_right=w-83;box_top=h-387;box_bottom=52
    c.setFillColor(HexColor('#dddddd'))
    c.rect(box_x,box_bottom,box_right-box_x,box_top-box_bottom,fill=1,stroke=0)
    c.setFillColor(HexColor('#111111'));c.setFont('AssignmentSerif',16.3)
    c.drawString(box_x+17,box_top-25,'9. Сфера как локально симметрическое пространство')
    c.setFont('AssignmentSerif',16.1)
    # First line indented as in the photo; S^2 has a true superscript.
    x=box_x+32;y=box_top-80
    prefix='Доказать, что сфера '
    c.drawString(x,y,prefix);x+=pdfmetrics.stringWidth(prefix,'AssignmentSerif',16.1)
    c.drawString(x,y,'S');x+=pdfmetrics.stringWidth('S','AssignmentSerif',16.1)
    c.setFont('AssignmentSerif',10);c.drawString(x,y+7,'2');x+=8
    c.setFont('AssignmentSerif',16.1)
    c.drawString(x,y,' с обычной метрикой — локально')
    c.drawString(box_x+5,y-25,'симметрическое пространство. Вычислить скалярную и гауссову')
    c.drawString(box_x+5,y-50,'кривизны')
    c.showPage();c.save()
    print(f'Built: {destination}')


if __name__=='__main__':main()
