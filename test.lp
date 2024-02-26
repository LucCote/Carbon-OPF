\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
 
Subject To
 R0: - C0 - C6 - C9 - C12 = -1
 R1: - C1 - C7 - C10 - C13 = 0
 R2: - C2 - C8 - C11 - C14 = -3
 R3: C6 <= 0
 R4: C7 <= 10
 R5: C8 <= 0
 R6: C9 <= 10
 R7: C10 <= 0
 R8: C11 <= 10
 R9: C12 <= 0
 R10: C13 <= 10
 R11: C14 <= 0
 R12: C6 = 0
 R13: - C3 + C4 + C7 = 0
 R14: C8 = 0
 R15: C3 - C4 + C9 = 0
 R16: C10 = 0
 R17: - C4 + C5 + C11 = 0
 R18: C12 = 0
 R19: C4 - C5 + C13 = 0
 R20: C14 = 0
 R21: C3 = 0
 R22: - C15 + C24 + C27 + C30 = 0
 R23: - C19 + C25 + C28 + C31 = 0
 R24: - C23 + C26 + C29 + C32 = 0
 R25: C16 = 0
 R26: C17 = 0
 R27: C20 = 0
 R28: C18 = 0
 R29: C21 = 0
 R30: C22 = 0
 R31: 1.25 C15 + 2 C16 + 2 C17 - 1.25 C24 - 2 C25 - 2 C26 >= 1
 R32: 1.25 C18 + 2 C19 + 2 C20 - 1.25 C27 - 2 C28 - 2 C29 >= 0
 R33: 1.25 C21 + 2 C22 + 2 C23 - 1.25 C30 - 2 C31 - 2 C32 >= 2
Bounds
 C0 <= 2
 C1 = 0
 C2 <= 2
 C3 >= -1000
 C4 >= -1000
 C5 >= -1000
 -1000 <= C6 <= 0
 -1000 <= C7 <= 10
 -1000 <= C8 <= 0
 -1000 <= C9 <= 10
 -1000 <= C10 <= 0
 -1000 <= C11 <= 10
 -1000 <= C12 <= 0
 -1000 <= C13 <= 10
 -1000 <= C14 <= 0
General Constraints
 GC0: C24 = MAX ( C6 , 0 )
 GC1: C25 = MAX ( C7 , 0 )
 GC2: C26 = MAX ( C8 , 0 )
 GC3: C27 = MAX ( C9 , 0 )
 GC4: C28 = MAX ( C10 , 0 )
 GC5: C29 = MAX ( C11 , 0 )
 GC6: C30 = MAX ( C12 , 0 )
 GC7: C31 = MAX ( C13 , 0 )
 GC8: C32 = MAX ( C14 , 0 )
End
