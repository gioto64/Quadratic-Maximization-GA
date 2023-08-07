import random
import matplotlib.pyplot as plot
from matplotlib.animation import FuncAnimation

class PlotData:
    def __init__(self):
        self.steps = []
        self.max = []
        self.average = []

        plot.ion()
        self.figure, self.ax = plot.subplots()
        self.line_max, = self.ax.plot(self.steps, self.max, label='Max Fitness')
        self.line_average, = self.ax.plot(self.steps, self.average, label='Average Fitness')
        self.ax.legend(loc='upper left')

    def append(self, i, max_fitness, average_fitness):
        self.steps.append(i)
        self.max.append(max_fitness)
        self.average.append(average_fitness)

        self.line_max.set_xdata(self.steps)
        self.line_max.set_ydata(self.max)
        self.line_average.set_xdata(self.steps)
        self.line_average.set_ydata(self.average)
        self.ax.relim()

        self.ax.autoscale_view()
        plot.draw()
        plot.pause(0.001)

    def final_show(self):
        plot.close()
        plot.plot(self.steps, self.max, label='Max Fitness')
        plot.plot(self.steps, self.average, label='Average Fitness')
        plot.legend(loc='upper left')
        plot.show(block=True)

        
class Cromozome:
    def __init__(self, shift, bits, len, precision, value):
        self.shift = shift
        self.bits = bits
        self.len = len
        self.precision = precision
        self.genes = [0] * bits
        self.from_value(value)

    def from_value(self, value):
        p2 = 2 ** (self.bits - 1)
        for i in range(self.bits - 1, -1, -1):
            if value >= p2:
                value -= p2
                self.genes[i] = 1
            p2 >>= 1
        return

    def to_value(self):
        value = 0
        p2 = 1
        for i in range(0, self.bits):
            value += self.genes[i] * p2
            p2 *= 2

        return value / (2 ** self.bits - 1) * self.len + self.shift

    def eval_fitness(self, function):
        value = self.to_value()
        return value * value * function[0] + value * function[1] + function[2]

    def print(self, file, function):
        print(f"Genes = {self.genes}", file = file)
        print(f"Value = {self.to_value()}", file = file)
        print(f"Fitness = {self.eval_fitness(function)}", file = file)

class CromozomeFactory:
    def __init__(self, l, r, precision):
        self.shift = l
        self.bits = 0;
        self.precision = 10 ** precision
        self.len = (r - l)

        p2 = 1
        tmp = self.len * self.precision 
        while p2 <= tmp:
            p2 = p2 * 2
            self.bits += 1
        
        return

    def create_cromozome(self, value):
        # value [l, r]
        value = (value - self.shift) #value [0, (r - l)]
        value = (value / self.len) * ((2 ** self.bits) - 1) #value [0, 2^bits - 1]
        cromozome = Cromozome(self.shift, self.bits, self.len, self.precision, value)
        return cromozome


def initialize_population(D, l, r, precision):
    cromozome_factory = CromozomeFactory(l, r, precision)
    population = []
    for _ in range(0, D):
        value = random.randint(l * (10 ** precision), r * (10 ** precision))
        value = value / (10 ** precision)
        population.append(cromozome_factory.create_cromozome(value))
    return (population, cromozome_factory)

def crossover(i, j, population):
    cut_off = random.randint(1, population[i].bits - 2)
    for k in range(cut_off, population[i].bits):
        population[i].genes[k], population[j].genes[k] = population[j].genes[k], population[i].genes[k]

    return

file = open("results.txt", "w")
D = int(input("Please input the population size: "))

a = float(input("Please input the coefficient of the quadratic term (X^2): "))
b = float(input("Please input the coefficient of the linear term (X): "))
c = float(input("Please input the coefficient of X: "))
function = (a, b, c) # f(X) = a * X^2 + b * X + c

l = float(input("Please input the lower bound of the domain: "))
r = float(input("Please input the upper bound of the domain: "))

precision = int(input("Please input the precision amount: "))
prob_crossover = float(input("Please input the probability of crossover: "))
prob_mutation = float(input("Please input the probability of mutation: "))
steps = int(input("Please input the number of steps the algorithm should take: "))

(population, cromozome_factory) = initialize_population(D, l, r, precision)

plot_data = PlotData()

mx_value = 0
x_value = 0

for i in range(0, steps):
    print(f"Current Step = {i}\n", file = file)

    total_fitness = 0
    fitness_values = []
    selection_probabilities = []
    acc_probabilities = []

    for cromozome in population:
        fitness_value = cromozome.eval_fitness(function)
        total_fitness += fitness_value
        fitness_values.append(fitness_value)

    selection_sum = 0
    acc_probabilities.append(selection_sum)

    max_fitness = fitness_values[0] 
    best_cromozome = 0
    for j in range(0, D):
        cromozome = population[j]
        if i == 0:
            cromozome.print(file, function)

        if fitness_values[j] > max_fitness:
            max_fitness = fitness_values[j]
            best_cromozome = j

        selection_probability = fitness_values[j] / total_fitness
        selection_sum += selection_probability

        selection_probabilities.append(selection_probability)
        acc_probabilities.append(selection_sum)
        if i == 0:
            print(f"Selection probability = {selection_probability}", file = file)
            print("", file = file)
    
    print(f"The average fitness is {total_fitness / D}\n", file = file)
    print(f"Best performer is cromozome {best_cromozome} with fitness {max_fitness}\n", file = file)

    if max_fitness > mx_value:
        mx_value = max_fitness
        x_value = population[best_cromozome].to_value() 

    plot_data.append(i, max_fitness, total_fitness / D)

    if i == 0:
        print(f"Cumulative probabilities = {acc_probabilities}", file = file)
        print("", file = file)
        print("Selection process begins:", file = file)

    new_population = []
    chosen_cromzomes = []

    new_population.append(population[best_cromozome])
    chosen_cromzomes.append(best_cromozome)
    
    for j in range(1, D):
        u = random.uniform(0, 1)
        if i == 0:
            print("", file = file)
            print(f"Generated number at step #{j} = {u}", file = file)
        (st, dr) = (0, D - 1)
        while st <= dr:
            mid = (st + dr) // 2
            if acc_probabilities[mid] < u:
                st = mid + 1
            else:
                dr = mid - 1

        chosen_cromzomes.append(st - 1)
        new_population.append(cromozome_factory.create_cromozome(population[st - 1].to_value()))
        
        if i == 0:
            print(f"Chosen cromozome at step #{j} = {st}", file = file)

    if i == 0:
        print("", file = file)
        print(f"Chosen cromozomes from old population = {chosen_cromzomes}", file = file)
        print(f"\nCrossover stage begins:\n", file = file)

    crossover_cromozomes = []
    # We start the range from one as we don't want to change the best performer
    for j in range(1, D):
        u = random.uniform(0, 1)
        if u < prob_crossover:
            crossover_cromozomes.append(j)

    if i == 0:
        print(f"Selected cromozomes for crossover = {crossover_cromozomes}", file = file)

    while len(crossover_cromozomes) > 1:
        j = random.randint(0, len(crossover_cromozomes) - 1)
        first_cromozome = crossover_cromozomes[j]
        crossover_cromozomes.remove(first_cromozome)

        j = random.randint(0, len(crossover_cromozomes) - 1)
        second_cromozome = crossover_cromozomes[j]
        crossover_cromozomes.remove(second_cromozome)

        if i == 0:
            print(f"Selected pair for crossover: {first_cromozome}, {second_cromozome}", file = file)
        
        crossover(first_cromozome, second_cromozome, new_population)

    if i == 0:
        print("\nPopulation after crossover:\n", file = file)
        for cromozome in new_population: 
            cromozome.print(file, function)
            print("", file = file)
        
    if i == 0:
        print("", file = file)
        print("Mutation process begins:\n", file = file)

    # We start the range from one as we don't want to change the best performer
    for k in range(1, D):
        cromozome = new_population[k]
        for j in range(0, cromozome.bits):
            u = random.uniform(0, 1)
            if u < prob_mutation:
                cromozome.genes[j] = cromozome.genes[j] ^ 1
                if i == 0:
                    print(f"Mutate the gene {j} of cromozome {k + 1} ", file = file)
    
    if i == 0:
        print("\nPopulation after mutations:\n", file = file)
        for cromozome in new_population: 
            cromozome.print(file, function)
            print("", file = file)

    population = new_population
    print(80 * '-', file = file)

plot_data.final_show()

print(f"The maximum value found is: f[{x_value}] = {mx_value}")

