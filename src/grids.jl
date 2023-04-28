export OneDGrid

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct OneDGrid

    dev :: AbstractDevice
    len  :: Int
    start :: Real
    stop  :: Real
    points :: LinRange
    step :: Real

    function OneDGrid( dev, len, start, stop )

        points = LinRange( start, stop, len+1 )[1:end-1]
        step = ( stop - start ) / len
        new( dev, len, start, stop, points, step)

    end

end
