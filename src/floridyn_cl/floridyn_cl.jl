# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function angSOWFA2world(deg_SOWFA)
    # Angle conversion SOWFA to world coordinates
    # deg_F = -deg_S + 270 deg
    # yaw angle defined clockwise, but for calculations counterclockwise
    deg_World = 270 - deg_SOWFA
    rad_World = deg2rad(deg_World)
    return rad_World
end
