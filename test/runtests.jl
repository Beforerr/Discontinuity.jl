using TestItems, TestItemRunner
@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(Discontinuity)
end

@testitem "Utility functions" begin
    @test Discontinuity.relative_variation([3 4 0; 4 3 0]) == 0.2
end
