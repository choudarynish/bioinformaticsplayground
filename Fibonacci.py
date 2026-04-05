# Function to get the nth Fibonacci number
def fibonacci_nth(n):
    if n <= 0:
        return "Invalid input: n must be a positive integer"
    elif n == 1:
        return 0  # First term
    elif n == 2:
        return 1  # Second term

    a, b = 0, 1
    for _ in range(3, n + 1):
        a, b = b, a + b
    return b

# Example usage
n = int(input("Enter n: "))  # Change this to the position you want
print(f"The {n}th Fibonacci number is:", fibonacci_nth(n))
