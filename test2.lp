\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  10 G1 + 15 G2 + 12 G3
Subject To
 R0: - G1 - F11 - F21 - F31 = -1
 R1: - G2 - F12 - F22 - F32 = 0
 R2: - G3 - F13 - F23 - F33 = -3
 R3: F11 <= 0
 R4: F12 <= 10
 R5: F13 <= 0
 R6: F21 <= 10
 R7: F22 <= 0
 R8: F23 <= 10
 R9: F31 <= 0
 R10: F32 <= 10
 R11: F33 <= 0
 R12: F11 = 0
 R13: - v1 + v2 + F12 = 0
 R14: F13 = 0
 R15: v1 - v2 + F21 = 0
 R16: F22 = 0
 R17: - v2 + v3 + F23 = 0
 R18: F31 = 0
 R19: v2 - v3 + F32 = 0
 R20: F33 = 0
 R21: v1 = 0
Bounds
 G1 <= 2
 G2 = 0
 G3 <= 2
 v1 >= -1000
 v2 >= -1000
 v3 >= -1000
 -1000 <= F11 <= 0
 -1000 <= F12 <= 10
 -1000 <= F13 <= 0
 -1000 <= F21 <= 10
 -1000 <= F22 <= 0
 -1000 <= F23 <= 10
 -1000 <= F31 <= 0
 -1000 <= F32 <= 10
 -1000 <= F33 <= 0
End
