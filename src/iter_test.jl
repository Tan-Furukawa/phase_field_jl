struct Fibo
end

# Base.iterate(::Fibo, state = (0, "hello")) =
# if state[1] < 0 && state[1] > 10
#   nothing
# else
#   (state[1] + 1, "hey")
:a
# end

Base.iterate(::Fibo, state=(0, 1)) =
       if state[1] < 0
       nothing
       else
       (state[1], (state[2], state[1] + state[2]))
       end

for x in Fibo()
  print(x, " ")
end

x = [1,2,3]
iterate(x,2)
iterate(x)

struct MyIterator
  state:: Tuple{Int, String}
  fin:: Int
end

function Base.iterate(itr::MyIterator, state = itr.state)
  if state[1] < 0 || state[1] > itr.fin
      return nothing
  end
  print(state)
  return (("a", state[1]^2), (state[1] + 1, "hello"))
end

itr = MyIterator((0, "10"), 3)
for s in itr
    println(s)
end
