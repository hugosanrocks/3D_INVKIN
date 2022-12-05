  #title
    ~/Desktop/cwp/bin/pslabel \
    t="Structured grid" size=50 \
    bcolor=white tcolor=black > labelreg.ps

    ~/Desktop/cwp/bin/pslabel \
    t="Non-structured grid" size=50 \
    bcolor=white tcolor=black > labelirreg.ps

    ~/Desktop/cwp/bin/psmerge \
    in=labelreg.ps translate=4,0.5 \
    in=regular.ps translate=0,0 \
    in=labelirreg.ps translate=19,0.5 \
    in=irregular.ps translate=8.37,-0.36 > out.ps

ps2pdf -dEPSCrop out.ps
pdfcrop --margins '-580 0 -650 0' out.pdf 
