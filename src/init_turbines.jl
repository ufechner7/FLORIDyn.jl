function turbineArrayProperties()
    # Foot position of the turbine x, y, z (NOT nacelle position, ground level!)
    Pos = [
        600   2400   0;  # T0
        1500  2400   0;  # T1
        2400  2400   0;  # T2
        600   1500   0;  # T3
        1500  1500   0;  # T4
        2400  1500   0;  # T5
        600   600    0;  # T6
        1500  600    0;  # T7
        2400  600    0;  # T8
    ]

    # Turbine type, all using "DTU 10MW"
    Type = [
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
        "DTU 10MW",
    ]

    # Initial states: a, yaw (deg), TI
    Init_States = [
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
        0.33  0  0.06;
    ]

    # Return as named tuple
    return (
        Pos = Pos,
        Type = Type,
        Init_States = Init_States
    )
end
