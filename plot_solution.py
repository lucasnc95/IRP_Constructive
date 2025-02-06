import sys
import os
import matplotlib.pyplot as plt
import networkx as nx

def plot_routes(filename):
    depots = {}
    customers = {}
    inventory_before = {}
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
        elif line.startswith("Inventory levels before delivery:"):
            section = "inventory_before"
        elif section == "depots":
            depot_info = list(map(float, line.strip().split()))
            depots[int(depot_info[0])] = (depot_info[1], depot_info[2])
        elif section == "customers":
            customer_info = list(map(float, line.strip().split()))
            customers[int(customer_info[0])] = (customer_info[1], customer_info[2])
        elif section == "inventory_before":
            customer_id, inventory = map(int, line.strip().split()[1:])
            inventory_before[customer_id] = inventory
        elif section == "period":
            if line.startswith("  Vehicle"):
                vehicle_id = int(line.split()[1][:-1])
                route = list(map(int, line.split(":")[1].strip().split()))
                routes[-1].append((vehicle_id, route))

    # Cores para os veículos (até 200 veículos)
    colors = plt.cm.tab20(range(200))

    # Pasta para salvar os plots
    plot_dir = 'Plots'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Função para determinar a cor do cliente baseado no nível de inventário
    def get_customer_color(inventory):
        if inventory < 0:
            return 'red'
        elif inventory > 0:
            return 'green'
        else:
            return 'yellow'

    # Plotagem dos gráficos de rotas para cada período
    for period in range(len(routes)):
        if not routes[period]:
            continue  # Pula períodos sem rotas

        G = nx.DiGraph()

        # Adiciona nós para depósitos e clientes
        for depot_id, (x, y) in depots.items():
            G.add_node(depot_id, pos=(x, y), color='black', shape='s')  # Depósito como quadrado preto
        for customer_id, (x, y) in customers.items():
            inventory = inventory_before.get(customer_id, 0)
            G.add_node(customer_id, pos=(x, y), color=get_customer_color(inventory), shape='o')

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
        node_shapes = [data['shape'] for _, data in G.nodes(data=True)]
        edge_colors = [data['color'] for _, _, data in G.edges(data=True)]
        edge_styles = [data['style'] if 'style' in data else 'solid' for _, _, data in G.edges(data=True)]
        edge_weights = [data['weight'] for _, _, data in G.edges(data=True)]

        plt.figure(figsize=(25, 25))

        # Desenho dos nós com formas específicas
        for shape in set(node_shapes):
            nx.draw_networkx_nodes(G, pos, nodelist=[n for n, d in G.nodes(data=True) if d['shape'] == shape],
                                   node_color=[d['color'] for n, d in G.nodes(data=True) if d['shape'] == shape],
                                   node_shape=shape, node_size=250, label=shape)
        
        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_weights, style=edge_styles)
        nx.draw_networkx_labels(G, pos, labels={n: str(n) for n in G.nodes()}, font_color='black')

        plt.title(f'Period {period} - Vehicle Routes')
        plt.legend(loc='best', markerscale=2, title="Node Shapes")
        plt.savefig(os.path.join(plot_dir, f'routes_period_{period}.png'))
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_filename")
    else:
        filename = sys.argv[1]
        plot_routes(filename)
