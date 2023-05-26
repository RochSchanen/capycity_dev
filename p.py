A = {"A0" : (841, 1189)}
for i in range(10):
    W, H = A[f"A{i}"]
    A[f"A{i+1}"] = (H/2, W)

print(A)