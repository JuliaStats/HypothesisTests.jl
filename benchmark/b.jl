using HypothesisTests, BenchmarkLite

type BTestProc <: Proc end
Base.string(::BTestProc) = "b-test"
Base.length(::BTestProc, cfg) = prod(cfg)
Base.isvalid(::BTestProc, cfg) = true
Base.start(::BTestProc, cfg) = (randn(cfg[1], cfg[2]), randn(cfg[1], cfg[2]))
Base.run(::BTestProc, cfg, s) = BTest(s[1], s[2])
Base.done(::BTestProc, cfg, s) = nothing

sizes = 4.^(2:5)
show(run(Proc[BTestProc()], vec([(x, y) for y in sizes, x in sizes])), unit=:msec)
