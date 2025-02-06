import os
import subprocess
import itertools

# Parâmetros das instâncias
a_values = [10, 25, 50, 100, 200]
b_values = [8, 18, 38]
c_values = [6, 9, 12]
prefixes = ["R_C", "R_RC", "R_R", "U_C", "U_RC", "U_R"]

# Número de soluções geradas por execução
n = 100  # valor de `n` a ser passado para o executável

# Número de vezes que cada instância será executada
executions_per_instance = 5  # Ajuste conforme necessário

# Valor de maxNoImprovement
max_no_improvement = 9  # Ajuste conforme necessário

# Caminho do executável e pasta das instâncias
executable = "./IRP"
instances_folder = "Instances"


# Gera todas as combinações possíveis de instâncias
for prefix, a, b, c in itertools.product(prefixes, a_values, b_values, c_values):
    instance_name = f"{prefix}_N{a}_Q{b}_T{c}"
    instance_path = os.path.join(instances_folder, f"{instance_name}.txt")
    
    # Verifica se a instância existe
    if not os.path.isfile(instance_path):
        print(f"Instância não encontrada: {instance_path}")
        continue
    
    # Executa a instância `executions_per_instance` vezes com sementes diferentes
    for run in range(executions_per_instance):
        seed = run + 1  # Define uma semente diferente para cada execução

        # Comando de execução
        command = [
            executable, 
            instance_path, 
            str(n),            # Número de soluções geradas por execução
            str(seed),         # Semente de randomização
            str(max_no_improvement)
        ]

        print(f"Executando: {' '.join(command)}")

        try:
            # Executa o comando e espera finalizar
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar {command}: {e}")
