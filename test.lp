\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  10 C0 + 15 C1 + 12 C2
Subject To
 R0: - C0 - C7 - C8 + C9 + C12 = -1
 R1: - C1 + C7 - C9 - C11 + C13 = 0
 R2: - C2 + C8 + C11 - C12 - C13 = -3
 R3: C6 <= 0
 R4: C7 <= 10
 R5: C8 <= 10
 R6: C9 <= 10
 R7: C10 <= 0
 R8: C11 <= 10
 R9: C12 <= 10
 R10: C13 <= 10
 R11: C14 <= 0
 R12: C6 = 0
 R13: - 1.428571428571429 C3 + 1.428571428571429 C4 + C7 = 0
 R14: - C3 + C5 + C8 = 0
 R15: 1.428571428571429 C3 - 1.428571428571429 C4 + C9 = 0
 R16: C10 = 0
 R17: - C4 + C5 + C11 = 0
 R18: C3 - C5 + C12 = 0
 R19: 1.428571428571429 C4 - 1.428571428571429 C5 + C13 = 0
 R20: C14 = 0
 R21: C3 = 0
Bounds
 C0 <= 2
 C1 <= 2
 C2 <= 3
 C6 = 0
 C7 <= 10
 C8 <= 10
 C9 <= 10
 C10 = 0
 C11 <= 10
 C12 <= 10
 C13 <= 10
 C14 = 0
End
