using Plots # unlike R but like Python: need to import most capabilities explicitly
#using LaTeXStrings # Plots can deal with LaTeX, but needs lots of escape characters; LaTeXStrings makes it easier
using Random

function double_well(x, r1, r2, r3, dt) # Like R, functions can return as if in the terminal
    x .- (x .- r1).*(x .- r2).*(x .- r3).*dt # dot notation is for vectorized operations
end

r1 = 1
r2 = 4
r3 = 7
dt = .01

T = 100

initialX = [10, 4.1, 3.9, 0.1]

X = zeros(Float64, T, length(initialX)) # a matrix of zeros with T rows and |initialX| columns
X[1, : ] = initialX # need the `:` operator to grab all columns, but otherwise slicing works as normal

for t in 2:T # sequences generated normally
    X[t, : ] = double_well(X[t-1, :], r1, r2, r3, dt)# .+ randn(length(initialX)).*0.1
end

plot(# https://docs.juliaplots.org/stable/
    1:T, X, lw = 2, 
    label = ["x_1" "x_2" "x_3" "x_4"],
    title = "Double Well System: Four Initial Values",
    xlabel = "t", ylabel = "X"
)
