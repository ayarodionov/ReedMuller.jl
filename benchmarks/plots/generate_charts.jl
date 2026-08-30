# Generates static SVG line charts (WER vs Eb/N0) for the RM(2,8-10)
# comparison in ../RESULTS.md, from the results hardcoded below (see
# benchmarks/results/compare_large_m8-10.log and
# benchmarks/results/compare_large_resume_m10.log for the full run
# that produced them — those logs are not tracked in the repo).
#
# Hand-rolled SVG (no plotting package dependency) using the
# dataviz-skill's 8-slot categorical palette and mark specs: 2px
# lines, filled circle markers, recessive gridlines, direct legend.
# Points with exactly 0 observed errors (out of 100,000 trials) are
# drawn as hollow markers at a 1e-5 floor with a footnote, since a
# true zero cannot be placed on a log axis.

using Printf

const EBN0 = [0.0, 1.0, 2.0, 3.0]
const FLOOR = 1e-5   # 1 / trials

# name => (color, [wer at 0,1,2,3 dB per code])
const SERIES = [
    ("Reed (hard)",           "#2a78d6"),
    ("Sidelnikov-Pershakov",  "#eb6834"),
    ("RPA",                   "#1baf7a"),
    ("Dumer FHT leaves",      "#eda100"),
    ("AED-8 / Dumer-FHT",     "#e87ba4"),
    ("DS16 FHT leaves",       "#008300"),
    ("GLP cyclic",            "#4a3aa7"),
    ("Graph search N=32",     "#e34948"),
]

const DATA = Dict(
    8 => Dict(
        "Reed (hard)"          => [1.00e+00, 9.99e-01, 9.95e-01, 9.69e-01],
        "Sidelnikov-Pershakov" => [5.38e-01, 2.02e-01, 3.18e-02, 1.21e-03],
        "RPA"                  => [1.81e-01, 3.56e-02, 2.62e-03, 7.00e-05],
        "Dumer FHT leaves"     => [6.73e-01, 4.15e-01, 1.76e-01, 4.13e-02],
        "AED-8 / Dumer-FHT"    => [2.63e-01, 6.93e-02, 7.76e-03, 2.50e-04],
        "DS16 FHT leaves"      => [2.75e-01, 8.41e-02, 1.34e-02, 8.60e-04],
        "GLP cyclic"           => [1.51e-01, 2.64e-02, 1.70e-03, 2.00e-05],
        "Graph search N=32"    => [1.20e-01, 1.89e-02, 1.10e-03, 1.00e-05],
    ),
    9 => Dict(
        "Reed (hard)"          => [1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00],
        "Sidelnikov-Pershakov" => [6.98e-01, 2.94e-01, 4.28e-02, 1.13e-03],
        "RPA"                  => [1.40e-01, 1.69e-02, 6.00e-04, 0.00e+00],
        "Dumer FHT leaves"     => [8.20e-01, 5.95e-01, 3.05e-01, 9.25e-02],
        "AED-8 / Dumer-FHT"    => [4.14e-01, 1.29e-01, 1.55e-02, 6.50e-04],
        "DS16 FHT leaves"      => [4.73e-01, 2.00e-01, 4.80e-02, 5.37e-03],
        "GLP cyclic"           => [2.08e-01, 3.75e-02, 2.15e-03, 2.00e-05],
        "Graph search N=32"    => [8.11e-02, 7.33e-03, 2.40e-04, 0.00e+00],
    ),
    10 => Dict(
        "Reed (hard)"          => [1.00e+00, 1.00e+00, 1.00e+00, 1.00e+00],
        "Sidelnikov-Pershakov" => [8.90e-01, 5.19e-01, 1.00e-01, 2.73e-03],
        "RPA"                  => [1.18e-01, 9.48e-03, 1.00e-04, 0.00e+00],
        "Dumer FHT leaves"     => [9.34e-01, 7.90e-01, 5.23e-01, 2.24e-01],
        "AED-8 / Dumer-FHT"    => [6.75e-01, 3.14e-01, 5.99e-02, 2.76e-03],
        "DS16 FHT leaves"      => [7.27e-01, 4.39e-01, 1.61e-01, 3.15e-02],
        "GLP cyclic"           => [4.21e-01, 1.16e-01, 1.00e-02, 1.10e-04],
        "Graph search N=32"    => [1.26e-01, 1.58e-02, 7.70e-04, 0.00e+00],
    ),
)

const W, H = 640, 420
const ML, MR, MT, MB = 64, 190, 24, 44   # margins: legend lives in the right margin
const PW, PH = W - ML - MR, H - MT - MB
const YMIN, YMAX = FLOOR, 1.0            # log10 axis range

xpix(ebn0) = ML + (ebn0 - EBN0[1]) / (EBN0[end] - EBN0[1]) * PW
ypix(wer) = MT + (log10(YMAX) - log10(max(wer, YMIN))) / (log10(YMAX) - log10(YMIN)) * PH

function chart_svg(m::Int)
    data = DATA[m]
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $W $H" font-family="-apple-system,Helvetica,Arial,sans-serif">""")
    println(io, """<rect width="$W" height="$H" fill="#fcfcfb"/>""")

    # Gridlines + y labels (log decades).
    for e in -5:0
        y = ypix(10.0^e)
        println(io, """<line x1="$ML" y1="$y" x2="$(ML+PW)" y2="$y" stroke="#dddad2" stroke-width="1"/>""")
        label = e == 0 ? "1" : "1e$e"
        println(io, """<text x="$(ML-8)" y="$(y+4)" font-size="12" fill="#52514e" text-anchor="end">$label</text>""")
    end
    # x gridlines + labels.
    for x in EBN0
        xp = xpix(x)
        println(io, """<line x1="$xp" y1="$MT" x2="$xp" y2="$(MT+PH)" stroke="#dddad2" stroke-width="1"/>""")
        println(io, """<text x="$xp" y="$(MT+PH+22)" font-size="12" fill="#52514e" text-anchor="middle">$(Int(x))</text>""")
    end
    println(io, """<text x="$(ML+PW/2)" y="$(H-8)" font-size="13" fill="#0b0b0b" text-anchor="middle">Eb/N0 (dB)</text>""")
    println(io, """<text x="16" y="$(MT+PH/2)" font-size="13" fill="#0b0b0b" text-anchor="middle" transform="rotate(-90 16 $(MT+PH/2))">Word error rate</text>""")
    println(io, """<text x="$ML" y="14" font-size="14" fill="#0b0b0b" font-weight="600">RM(2, $m) — n=$(2^m), BI-AWGN, 100,000 trials/point</text>""")

    # Series lines + markers.
    for (name, color) in SERIES
        vals = data[name]
        pts = [(xpix(EBN0[i]), ypix(vals[i])) for i in eachindex(EBN0)]
        path = join(["$(i==1 ? "M" : "L")$(round(x,digits=1)),$(round(y,digits=1))" for (i,(x,y)) in enumerate(pts)], " ")
        println(io, """<path d="$path" fill="none" stroke="$color" stroke-width="2.5"/>""")
        for (i, (x, y)) in enumerate(pts)
            if vals[i] <= 0.0
                println(io, """<circle cx="$x" cy="$y" r="5" fill="#fcfcfb" stroke="$color" stroke-width="2.5"/>""")
            else
                println(io, """<circle cx="$x" cy="$y" r="4.5" fill="$color"/>""")
            end
        end
    end

    # Legend (right margin).
    lx = ML + PW + 16
    ly = MT + 8
    for (name, color) in SERIES
        println(io, """<circle cx="$(lx+5)" cy="$(ly+4)" r="5" fill="$color"/>""")
        println(io, """<text x="$(lx+16)" y="$(ly+8)" font-size="12" fill="#0b0b0b">$name</text>""")
        ly += 20
    end
    println(io, """<text x="$lx" y="$(ly+8)" font-size="10.5" fill="#52514e">○ = 0 errors observed</text>""")
    println(io, """<text x="$lx" y="$(ly+22)" font-size="10.5" fill="#52514e">(floor shown at 1/trials)</text>""")

    println(io, "</svg>")
    String(take!(io))
end

for m in (8, 9, 10)
    open(joinpath(@__DIR__, "rm2_$(m).svg"), "w") do f
        write(f, chart_svg(m))
    end
    println("wrote rm2_$(m).svg")
end
