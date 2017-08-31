(TeX-add-style-hook
 "lab4report"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "twocolumn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "inputenc"
    "amsmath"
    "amsthm"
    "amssymb"
    "subcaption"
    "graphicx"
    "listings"
    "color")
   (TeX-add-symbols
    "listingsfont"
    "listingsfontinline")
   (LaTeX-add-labels
    "sec:code-overview"
    "sec:disc-part"
    "sec:execution"
    "sec:results"
    "sec:code")))

