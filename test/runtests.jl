tests = [
    f for f in readdir(@__DIR__) if !(isnothing(match(r"test_", f))) && !(isnothing(match(r".jl", f)))
]

println("Starting tests")
for t in ["test_days_ordering.jl"]
  println("\n","#"^80)
  println("Running ", t,)
  println("#"^80, "\n")
  include(t)
end
