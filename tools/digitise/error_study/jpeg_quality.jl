# Estimate the JPEG quality of the PMC figure images from their quantisation tables.
# Reads the DQT and SOF markers, and finds the IJG (libjpeg) quality q whose scaled standard
# luminance / chrominance tables match the file's tables best. Also reports chroma subsampling.
#   julia --project=tools/digitise tools/digitise/error_study/jpeg_quality.jl <jpg>...
using Printf

const STD_LUM = [16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 14 13 16 24 40 57 69 56; 14 17 22 29 51 87 80 62;
                 18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 49 64 78 87 103 121 120 101; 72 92 95 98 112 100 103 99]
const STD_CHR = [17 18 24 47 99 99 99 99; 18 21 26 66 99 99 99 99; 24 26 56 99 99 99 99 99; 47 66 99 99 99 99 99 99;
                 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99; 99 99 99 99 99 99 99 99]
# zigzag order -> natural (row-major) index
const ZZ = [0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18, 11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28,
            35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63]

ijg(T, q) = (s = q < 50 ? 5000 ÷ q : 200 - 2q; clamp.((T .* s .+ 50) .÷ 100, 1, 255))

function jpeg_info(path)
    b = read(path)
    tables = Dict{Int, Matrix{Int}}()
    sampling = Tuple{Int, Int, Int}[]
    H = W = 0
    i = 3
    while i < length(b)
        b[i] == 0xff || (i += 1; continue)
        m = b[i + 1]
        (m == 0xd8 || m == 0x01 || 0xd0 <= m <= 0xd7) && (i += 2; continue)
        m == 0xda && break
        L = Int(b[i + 2]) << 8 | b[i + 3]
        seg = b[i + 4:i + 1 + L]
        if m == 0xdb
            k = 1
            while k <= length(seg)
                pq, tq = seg[k] >> 4, seg[k] & 0x0f
                vals = pq == 0 ? Int.(seg[k + 1:k + 64]) : [Int(seg[k + 2j - 1]) << 8 | seg[k + 2j] for j in 1:64]
                Q = zeros(Int, 8, 8)
                for (z, v) in zip(ZZ, vals); Q[z ÷ 8 + 1, z % 8 + 1] = v; end
                tables[tq] = Q
                k += 1 + 64 * (pq + 1)
            end
        elseif m in (0xc0, 0xc1, 0xc2)
            H = Int(seg[2]) << 8 | seg[3]; W = Int(seg[4]) << 8 | seg[5]
            nc = seg[6]
            for c in 1:nc
                s = seg[6 + 3(c - 1) + 2]
                push!(sampling, (Int(seg[6 + 3(c - 1) + 1]), Int(s >> 4), Int(s & 0x0f)))
            end
        end
        i += 2 + L
    end
    fit(T, Q) = (qs = 1:100; e = [sum(abs.(ijg(T, q) .- Q)) for q in qs]; (qs[argmin(e)], minimum(e)))
    ql, el = haskey(tables, 0) ? fit(STD_LUM, tables[0]) : (0, -1)
    qc, ec = haskey(tables, 1) ? fit(STD_CHR, tables[1]) : (0, -1)
    sub = length(sampling) >= 3 ? (sampling[1][2] == 2 && sampling[1][3] == 2 ? "4:2:0" : sampling[1][2] == 2 ? "4:2:2" : "4:4:4") : "grey"
    return (; H, W, q_lum = ql, err_lum = el, q_chr = qc, err_chr = ec, subsampling = sub)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("file,width,height,quality_luminance,table_mismatch_lum,quality_chrominance,table_mismatch_chr,subsampling")
    for p in ARGS
        r = jpeg_info(p)
        println(join((basename(p), r.W, r.H, r.q_lum, r.err_lum, r.q_chr, r.err_chr, r.subsampling), ","))
    end
end
