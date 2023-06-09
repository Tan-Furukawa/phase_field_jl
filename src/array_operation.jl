module ArrayOperation

export reject_element!
function reject_element!(arr::AbstractArray, tpl::Tuple)
  a, b = tpl
  arr[arr .< a] .= a
  arr[arr .> b] .= b
  return arr
end

end