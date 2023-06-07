function delete_png_files(folder_path)
  files = readdir(folder_path)
  for file in files
      if endswith(file, ".png")
          file_path = joinpath(folder_path, file)
          rm(file_path)
          println("Deleted file: $file_path")
      end
  end
end

# 削除するフォルダのパスを指定して実行する例
folder_path = "result/png"
delete_png_files(folder_path)