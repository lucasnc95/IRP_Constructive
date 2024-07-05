import sys
import matplotlib.pyplot as plt
import numpy as np

def plot_routes(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    periods = []
    routes = {}
    inventory_levels = {}
    period = -1

    for line in lines:
        if line.startswith("Period"):
            period = int(line.split()[1][:-1])
            periods.append(period)
            routes[period] = []
        elif line.startswith("  Vehicle"):
            vehicle_id = int(line.split()[1][:-1])
            route = []
            for pair in line.split(":")[1].strip().split(" "):
                customer_id, amount = pair.strip("()").split(",")
                route.append((int(customer_id), int(amount)))
            routes[period].append(route)
        elif line.startswith("Inventory levels"):
            break

    inventory_lines = lines[lines.index("Inventory levels:\n") + 1:]
    period = -1

    for line in inventory_lines:
        if line.startswith("Period"):
            period = int(line.split()[1][:-1])
            inventory_levels[period] = []
        else:
            customer_id = int(line.split()[1][:-1])
            inventory = int(line.split(":")[1].strip())
            inventory_levels[period].append((customer_id, inventory))

    # Plot routes
    plt.figure(figsize=(12, 8))
    for period in periods:
        plt.subplot(1, len(periods), period + 1)
        for vehicle_route in routes[period]:
            x_coords = []
            y_coords = []
            for customer_id, _ in vehicle_route:
                if customer_id == 0:
                    x, y = 0, 0  # Coordinates for the depot
                else:
                    customer = next(c for c in irp.customers if c.id == customer_id)
                    x, y = customer.x, customer.y
                x_coords.append(x)
                y_coords.append(y)
            plt.plot(x_coords, y_coords, marker='o')
        plt.title(f'Period {period}')
        plt.xlabel('X')
        plt.ylabel('Y')

    plt.tight_layout()
    plt.savefig('routes.png')

    # Plot inventory levels
    plt.figure(figsize=(12, 8))
    for period in periods:
        customer_ids = [c[0] for c in inventory_levels[period]]
        inventories = [c[1] for c in inventory_levels[period]]
        plt.bar(customer_ids, inventories, label=f'Period {period}')

    plt.xlabel('Customer ID')
    plt.ylabel('Inventory Level')
    plt.legend()
    plt.title('Inventory Levels Over Time')
    plt.savefig('inventory_levels.png')

if __name__ == "__main__":
    filename = sys.argv[1]
    plot_routes(filename)
