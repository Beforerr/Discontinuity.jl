"""the MVAB method can achieve acceptable accuracy when either |B|/|B| > 0.05 or ω > 60°. @liuFailuresMinimumVariance2023"""
assign_accuracy!(df) = @transform!(df, :accuracy = (:ω .> 60) .| (:dBmag_over_Bmag .> 0.05))
filter_low_accuracy(df) = filter(:accuracy => ==(true), df)