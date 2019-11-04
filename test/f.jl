using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "F-tests" begin

    @testset "Basic F-test" begin
        Random.seed!(123)
        y1_h0 = 4 .+ randn(500)
        y2_h0 = 4 .+ randn(400)

        t = FTest(y1_h0, y2_h0)

        @test t.n1 == 500
        @test t.n2 == 400
        @test t.F ≈ 1.234701 rtol = 1e-4
        @test pvalue(t) ≈ 0.027496 rtol = 1e-4
        @test default_tail(t) == :both

        y1_h1 = 0.8*randn(200)
        y2_h1 = 1.2*randn(120)

        t = FTest(y1_h1, y2_h1)

        @test t.n1 == 200
        @test t.n2 == 120
        @test t.F ≈ 0.436056 rtol = 1e-4
        @test pvalue(t) < 1e-5
        @test default_tail(t) == :both
        @test pvalue(t; tail = :left) < 1e-5
        @test pvalue(t; tail = :right) > 0.999
    end

    @testset "Homoscedasticity" begin
        y_h0 = [2.33982, 1.38455, -0.16677, -0.495884, -1.19301, -1.4365, 1.87554, -0.874521, 
            0.558398, 1.43168, -0.325775, -1.127, 0.381937, 0.309022, -1.48125, -0.848039, 
            0.107054, -0.141358, -0.973084, -0.424127, -2.0171, 0.659719, -0.605461, 0.60834, 
            0.301879, 1.31643, -0.111021, 1.12901, 0.179706, 0.558483, -1.4038, 1.46825, 
            -0.444613, 0.609621, 0.420681, 0.203042, 0.0243958, 1.31061, -0.624273, 0.190979, 
            0.0453361, -0.149526, -0.995021, 1.21091, 0.248534, -0.359841, 0.955801, 0.746285, 
            0.673524, -1.42755, 0.905993, -0.55439, 1.46732, 1.37445, -0.0448962, 1.08624, 0.435137, 
            -1.06591, 0.101213, -0.136897, -2.469, -1.30757, 0.00990314, -0.939588, -1.31599, 0.410133, 
            -1.13285, -0.66851, 0.957774, 0.31636, 0.796974, -2.50988, 1.79729, -0.996714, -0.00277338, 
            0.889443, -0.256234, -0.280391, 0.330453, 0.525693, -0.995728, 0.355456, 0.331616, 0.152123, 
            0.0817606, -1.24144, -0.286701, 0.966459, -0.239041, -0.261518,-0.160022, -0.541692, 0.421619, 
            1.611, 0.573377, 0.182253, -2.40633, -1.83369, 0.786375, -1.14539, 0.991982, 0.471455, 1.64664, 
            -0.437902, -1.4579, 0.137847, 0.545115, 1.81469, -0.650005, -0.152165, 1.57597, 0.49686, 
            1.10803, -0.806275, -1.67866, -0.5121, 0.212479, -0.237583, 1.5334, 0.489391, -0.767837, 
            -0.541198, 0.60911, 2.36357, -0.496155, -2.36264, -1.06455, 0.0244807, -0.532804, -0.830876, 
            -0.66738, -1.38578, 0.945989, 0.262255, 0.415198, 1.1703, -0.489673, -0.528454, 2.06919, 
            -0.322042, -1.25894, 0.0338996, -0.616511, -2.29975, 0.268315, -0.881625, 0.114956, 0.00810337, 
            -0.989782, -0.810291, -0.470241, -0.587929, 1.18924, 0.233295, 0.379715, -0.268514, -0.517191, 
            -0.240446, -1.35477, -1.398, -0.380054, 2.02965, -0.39639, -0.635511, 1.92757, 0.0937969, 
            1.28987, 2.24724, 0.184783, 1.62242, -0.687497, 1.83239, -0.343518, -2.44318, -0.257667, 
            0.741863, 0.145903, -0.450616, 0.384291, 0.0129242, -1.11626, 1.45973, -1.27855, -0.154354, 
            1.67222, -0.685274, -0.742184, 0.488662, 0.805217, -1.22203, -1.11926, 1.29633, -0.415697, 
            -0.319425, -1.14136, -0.59144, -0.350336, -0.426482, 0.275284, 1.58645]

        t = HomoscedasticityFTest(y_h0, 100)

        @test t.n1 == 100
        @test t.n2 == 100
        @test t.F ≈ 1.11576 rtol = 1e-4
        @test pvalue(t) ≈ 0.58681 rtol = 1e-4
        @test default_tail(t) == :both

        y_h1 = [-1.01244, -0.551161, -0.373784, -0.907554, 0.341815, 1.2635, 0.439993, 0.191919, -0.0450509,
            -0.155678, -0.329643, -1.1579, -0.109133, 0.849377, -0.29017, -0.457472, -1.85472, -0.672236, 
            -1.1841, 0.769306, -1.0916, -0.869373, 1.1198, -0.428124, -0.499798, 0.366467, -0.567023, 
            -0.106883, 0.153275, 0.211471, -0.108504, -0.0261862, -0.120691, -0.206043, 0.669944, -0.933982, 
            0.422045, -0.582051, -0.354932, 0.3276, -0.20954, 1.02069, -0.828192, -0.404579,-0.538621, 0.175236, 
            -0.0875101, 1.82075, 1.51047, 0.356628, -1.86363, -0.268291, 0.811953, -0.222526, 1.5156, -1.0736, 
            -0.729931, -0.0572633, -0.16562, -0.567673, 0.608075, -0.163243, 0.385735, 1.02767, 0.396526, 
            0.0903397, 0.0159794, 0.523044, 0.373952, 0.501596, 0.654368, 0.703847, 1.41553, 0.0718486, 
            -0.375002, -1.26181, 0.991733, 0.358758, 0.311868, -0.855431, -1.3727, 0.750549, 0.313968, 0.272195, 
            -0.912945, 0.321617, 0.736034, 0.0882905, 0.387182, 1.17084, -1.3751, 0.0562779, -0.905189, -0.363472, 
            -1.13253, -0.828526, -1.4749, 0.490708, 0.548871, -0.182271, 1.26973, -0.109062, 1.23201, 1.46129, 
            -1.89326, -0.0225106, -0.38015, -0.321539, -1.54044, 1.22308, 1.54552, 1.25847, -0.361827, -3.24469, 
            0.648251, -0.377657, 2.51312, -0.263548, 1.40933,-1.16114, -0.247717, -2.01171, 1.35864, 0.903761, 
            0.397161, 0.872612, 0.987669, -0.979173, 0.795435, -1.98066, -0.316685, -0.808121, 0.773451, -2.29368, 
            -0.144789, -1.77911, 0.619739, -0.298868, 0.709764, -0.0311074, 0.11051, 0.508197, -1.6656, -0.196769, 
            0.113227, -0.041016, -0.0199208, 2.17162, 1.39875, -2.852, 1.31147, 0.509837, 0.715059, 1.68435, 
            1.27134, 1.44392, 0.0464478, 1.90758, -1.15398, -3.08975, -0.136464, 1.36812, -0.14568, -1.18536, 
            -0.5973, 1.44578, -1.08278, 1.50277, 2.17551, 0.187916, -0.288632, -1.09679, -0.588233, -1.19962, 
            0.987613, 0.252459, 0.273592, 1.3616, 0.719356, 0.403357, 1.13866, 2.04045, 1.58365, -0.480208, 
            1.90633, 0.193938, -1.11014, -2.12503,-1.77391, 0.478253, 0.804435, 1.33476, -1.42924, 0.697926, 
            -0.347439, 0.0500163, -0.41421, 0.972537, 1.66392, 0.818117]

        t = HomoscedasticityFTest(y_h1, 80)

        @test t.n1 == 80
        @test t.n2 == 80
        @test t.F ≈ 2.61080 rtol = 1e-4
        @test pvalue(t) ≈ 2.99304e-5 rtol = 1e-4
        @test default_tail(t) == :both

        @test_throws ArgumentError HomoscedasticityFTest(y_h1, 101)

        @test pvalue(t; tail = :left) > 0.999
        @test pvalue(t; tail = :right) < 2e-5
    end

end