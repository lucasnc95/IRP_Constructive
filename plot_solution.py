import sys
import os
import matplotlib.pyplot as plt
import networkx as nx

def plot_routes(filename):
    depots = {}
    customers = {}
    routes = []

    with open(filename, 'r') as file:
        lines = file.readlines()

    section = None
    for line in lines:
        if line.startswith("Depots:"):
            section = "depots"
        elif line.startswith("Customers:"):
            section = "customers"
        elif line.startswith("Period"):
            section = "period"
            period = int(line.split()[1][:-1])
            routes.append([])
        elif line.startswith("Inventory levels:"):
            section = "inventory"
        elif section == "depots":
            depot_info = list(map(float, line.strip().split()))
            depots[int(depot_info[0])] = (depot_info[1], depot_info[2])
        elif section == "customers":
            customer_info = list(map(float, line.strip().split()))
            customers[int(customer_info[0])] = (customer_info[1], customer_info[2])
        elif section == "period":
            if line.startswith("  Vehicle"):
                vehicle_id = int(line.split()[1][:-1])
                route = list(map(int, line.split(":")[1].strip().split()))
                routes[-1].append((vehicle_id, route))

    # Cores para os veículos (até 20 veículos)
    colors = plt.cm.tab20(range(20))

    # Pasta para salvar os plots
    plot_dir = 'Plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Plotagem dos gráficos de rotas para cada período
    for period in range(len(routes)):
        if not routes[period]:
            continue  # Pula períodos sem rotas

        G = nx.DiGraph()

        # Adiciona nós para depósitos e clientes
        for depot_id, (x, y) in depots.items():
            G.add_node(depot_id, pos=(x, y), color='red')
        for customer_id, (x, y) in customers.items():
            G.add_node(customer_id, pos=(x, y), color='blue')

        # Adiciona arestas para as rotas dos veículos
        for i, (vehicle_id, route) in enumerate(routes[period]):
            if len(route) > 2 or (len(route) == 2 and route[0] != route[1]):  # Ignora rotas 0 -> 0
                for j in range(len(route) - 1):
                    if route[j] == 0 or route[j + 1] == 0:  # Arestas que envolvem o depósito
                        G.add_edge(route[j], route[j + 1], color=colors[i % 20], weight=0.5, style='dotted')
                    else:
                        G.add_edge(route[j], route[j + 1], color=colors[i % 20], weight=2)

        pos = nx.get_node_attributes(G, 'pos')
        node_colors = [data['color'] for _, data in G.nodes(data=True)]
        edge_colors = [data['color'] for _, _, data in G.edges(data=True)]
        edge_styles = [data['style'] if 'style' in data else 'solid' for _, _, data in G.edges(data=True)]
        edge_weights = [data['weight'] for _, _, data in G.edges(data=True)]

        plt.figure(figsize=(25, 25))
        nx.draw(G, pos, node_color=node_colors, edge_color=edge_colors, width=edge_weights, with_labels=True,
                node_size=250, font_size=10, font_color='black', arrowsize=15, style=edge_styles,
                connectionstyle='arc3, rad=0.0')
        plt.title(f'Period {period} - Vehicle Routes')
        plt.savefig(os.path.join(plot_dir, f'routes_period_{period}.png'))
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_filename")
    else:
        filename = sys.argv[1]
        plot_routes(filename)
