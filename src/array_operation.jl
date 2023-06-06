module array_operation

export reject_element
function reject_element(arr::AbstractArray, tpl::Tuple)
  a, b = tpl
  new_arr = deepcopy(arr)
  new_arr[new_arr .< a] .= a
  new_arr[new_arr .> b] .= b
  return new_arr
end

end