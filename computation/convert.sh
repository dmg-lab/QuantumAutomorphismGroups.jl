##!/bin/bash
mutable struct ClosingOrder
    position::Position
    volume::Real
    fulfilled::Bool
    fulfillment_date::Date
end


# Define the filename and the pattern to search for
filename="test.tex"
pattern="\\begin{tabular}"

begin_tabular_line=$(sed -n '/\\begin{tabular}/p' "$filename")

# Print the captured line
echo "Captured line: $begin_tabular_line"


newtext="\\begin{tabular}, \\end{tabular}"

sed -i '/\\begin{tabular}/,/\\end{tabular}/c\\test' test.tex
cat test.tex


sed -i '/\\begin{tabular}{.*}\n.*?\\hline\n/,/\\label{Tab:computational-results}/c\\Dummytextreplacement' test.tex
