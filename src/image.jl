using Plots

function sort_file_names(names)
  arr = map(x -> parse(Int, match(r"\d+", x).match), names)
  sorted_names = names[sortperm(arr)]
  return sorted_names
end

anim = Animation()
png_files = readdir("result/png", join=true)
png_files = sort_file_names(png_files)

# PNGファイルを読み込み、アニメーションフレームとして格納
frames = [load(png_file) for png_file in png_files]
println("making plot")
for f in frames
  plt = plot(f)
  frame(anim, plt)
  print(".")
end

# GIFに変換して保存
# gif_path = "result.gif"
gif(anim, "result/gif/test.gif", fps = 30)