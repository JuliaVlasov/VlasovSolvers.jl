export OneDGrid

struct OneDGrid

    dev :: AbstractDevice
    len  :: Int
    start :: Real
    stop  :: Real
    points :: LinRange

    function OneDGrid( dev, len, start, stop )

        points = LinRange( start, stop, len )

        new( dev, len, start, stop, points)

    end

end
