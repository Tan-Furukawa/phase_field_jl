# struct Fibo
# end

# # Base.iterate(::Fibo, state = (0, "hello")) =
# # if state[1] < 0 && state[1] > 10
# #   nothing
# # else
# #   (state[1] + 1, "hey")
# # end

# Base.iterate(::Fibo, state=(0, 1)) =
#        if state[1] < 0
#        nothing
#        else
#        (state[1], (state[2], state[1] + state[2]))
#        end

# for x in Fibo()
#   print(x, " ")
# end

# x = [1,2,3]
# iterate(x,2)
# iterate(x)

# struct MyIterator
#   state:: Tuple{Int, String}
#   fin:: Int
# end

# function Base.iterate(itr::MyIterator, state = itr.state)
#   if state[1] < 0 || state[1] > itr.fin
#       return nothing
#   end
#   print(state)
#   return (("a", state[1]^2), (state[1] + 1, "hello"))
# end

# itr = MyIterator((0, "10"), 3)
# for s in itr
#     println(s)
# end


struct MyIterator
  data::Vector{Int}
  index::Int
end


function Base.iterate(iter::MyIterator, state=iter.index)
  if state > length(iter.data)
      return nothing
  end

  next_element = iter.data[state]
  next_state = state + 1

  return (next_element, next_state)
end

data = [1, 2, 3, 4, 5]
my_iter = MyIterator(data, 1)

for element in my_iter
  println(element)
end

