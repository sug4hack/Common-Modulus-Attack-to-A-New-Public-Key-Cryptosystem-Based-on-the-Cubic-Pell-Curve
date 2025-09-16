import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


plt.style.use('default')
sns.set_palette("Set2")

df = pd.read_csv('/home/dandi/smt8/papernew/eksperimen/cpc_attack_bench.csv')

def basic_analysis(df):
    print("=== BASIC DATASET INFO ===")
    print(f"Dataset size: {len(df)}")
    print(f"Tested modulus sizes: {sorted(df['N_bits'].unique())} bits")
    print(f"Trials per modulus size: {df.groupby('N_bits').size().tolist()}")
    
    # Calculate success rate
    success_rate = df.groupby('N_bits')['success'].mean() * 100
    print("\n=== SUCCESS RATE ===")
    for bits, rate in success_rate.items():
        print(f"N={bits} bits: {rate:.2f}% success")
    
    return success_rate

def runtime_analysis(df):
    print("\n=== ATTACK RUNTIME ANALYSIS ===")
    
    attack_cols = ['t_recover_a_s', 't_cma_s']
    attack_labels = ['Recover a', 'CMA Attack']
    
    runtime_stats = df.groupby('N_bits')[attack_cols].agg(['mean', 'std', 'min', 'max'])
    
    for bits in sorted(df['N_bits'].unique()):
        print(f"\n--- Runtime for N={bits} bits ---")
        bit_data = df[df['N_bits'] == bits]
        for i, col in enumerate(attack_cols):
            mean_val = bit_data[col].mean()
            std_val = bit_data[col].std()
            print(f"{attack_labels[i]}: {mean_val:.6f} ± {std_val:.6f} seconds")
        
        total_attack_time = bit_data[attack_cols].sum(axis=1)
        print(f"Total attack time: {total_attack_time.mean():.6f} ± {total_attack_time.std():.6f} seconds")
    
    return runtime_stats

def memory_analysis(df):
    print("\n=== MEMORY ANALYSIS ===")

    memory_stats = df.groupby('N_bits')['peak_mem_kb'].agg(['mean', 'std', 'min', 'max'])
    
    for bits in sorted(df['N_bits'].unique()):
        print(f"\n--- Memory usage for N={bits} bits ---")
        bit_data = df[df['N_bits'] == bits]
        mean_mem = bit_data['peak_mem_kb'].mean()
        std_mem = bit_data['peak_mem_kb'].std()
        mean_mem_mb = mean_mem / 1024
        print(f"Peak memory: {mean_mem:.2f} ± {std_mem:.2f} KB ({mean_mem_mb:.2f} MB)")
    
    return memory_stats

def complexity_analysis(df):
    print("\n=== ATTACK COMPLEXITY ANALYSIS ===")
    
    n_bits = sorted(df['N_bits'].unique())
    
    # Calculate total attack time (recover a + CMA)
    df['t_attack_s'] = df['t_recover_a_s'] + df['t_cma_s']
    
    # Average time and memory per modulus size
    avg_times = df.groupby('N_bits')['t_attack_s'].mean()
    avg_memory = df.groupby('N_bits')['peak_mem_kb'].mean()
    
    # Log-transform for growth analysis
    log_bits = np.log(n_bits)
    log_times = np.log(avg_times.values)
    log_memory = np.log(avg_memory.values)
    
    # Linear regression for attack time
    time_slope, time_intercept, time_r_value, time_p_value, time_std_err = \
        stats.linregress(log_bits, log_times)
    
    # Linear regression for memory
    mem_slope, mem_intercept, mem_r_value, mem_p_value, mem_std_err = \
        stats.linregress(log_bits, log_memory)
    
    print(f"Attack time complexity: O(n^{time_slope:.4f}) with R² = {time_r_value**2:.4f}")
    print(f"Memory complexity: O(n^{mem_slope:.4f}) with R² = {mem_r_value**2:.4f}")
    
    return {
        'time_complexity': time_slope,
        'memory_complexity': mem_slope,
        'time_r_squared': time_r_value**2,
        'memory_r_squared': mem_r_value**2
    }

def create_visualizations(df, complexity_results):
    n_bits = sorted(df['N_bits'].unique())
    
    plt.figure(figsize=(10, 6))
    success_rate = df.groupby('N_bits')['success'].mean() * 100
    plt.bar(n_bits, success_rate.values)
    plt.xlabel('Modulus Size (bits)')
    plt.ylabel('Success Rate (%)')
    plt.title('Attack Success Rate Based on Modulus Size')
    for i, v in enumerate(success_rate.values):
        plt.text(n_bits[i], v + 1, f'{v:.1f}%', ha='center')
    plt.savefig('success_rate.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    plt.figure(figsize=(12, 8))
    attack_cols = ['t_recover_a_s', 't_cma_s']
    attack_labels = ['Recover a', 'CMA Attack']
    colors = plt.cm.Set2(np.linspace(0, 1, len(attack_cols)))
    
    for i, col in enumerate(attack_cols):
        avg_times = df.groupby('N_bits')[col].mean()
        plt.plot(n_bits, avg_times.values, 'o-', label=attack_labels[i], color=colors[i], linewidth=2)
    
    df['t_attack_total'] = df['t_recover_a_s'] + df['t_cma_s']
    avg_total = df.groupby('N_bits')['t_attack_total'].mean()
    plt.plot(n_bits, avg_total.values, 's-', label='Total Attack', color='red', linewidth=2)
    
    plt.xlabel('Modulus Size (bits)')
    plt.ylabel('Time (seconds)')
    plt.title('Attack Operation Runtime Based on Modulus Size')
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig('attack_runtime.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    plt.figure(figsize=(10, 6))
    avg_memory = df.groupby('N_bits')['peak_mem_kb'].mean() / 1024  # Convert to MB
    plt.plot(n_bits, avg_memory.values, 's-', color='purple', linewidth=2)
    plt.xlabel('Modulus Size (bits)')
    plt.ylabel('Peak Memory (MB)')
    plt.title('Memory Usage Based on Modulus Size')
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig('memory_usage.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    log_bits = np.log(n_bits)
    df['t_attack_total'] = df['t_recover_a_s'] + df['t_cma_s']
    log_times = np.log(df.groupby('N_bits')['t_attack_total'].mean().values)
    ax1.scatter(log_bits, log_times, s=80)
    
    time_slope = complexity_results['time_complexity']
    time_intercept = np.mean(log_times) - time_slope * np.mean(log_bits)
    regression_line = time_slope * log_bits + time_intercept
    ax1.plot(log_bits, regression_line, 'r-', 
             label=f'Slope: {time_slope:.3f}\nR²: {complexity_results["time_r_squared"]:.3f}')
    
    ax1.set_xlabel('log(Modulus Size)')
    ax1.set_ylabel('log(Attack Time)')
    ax1.set_title('Attack Time Complexity Analysis')
    ax1.legend()
    
    log_memory = np.log(df.groupby('N_bits')['peak_mem_kb'].mean().values)
    ax2.scatter(log_bits, log_memory, s=80)

    mem_slope = complexity_results['memory_complexity']
    mem_intercept = np.mean(log_memory) - mem_slope * np.mean(log_bits)
    regression_line = mem_slope * log_bits + mem_intercept
    ax2.plot(log_bits, regression_line, 'r-', 
             label=f'Slope: {mem_slope:.3f}\nR²: {complexity_results["memory_r_squared"]:.3f}')
    
    ax2.set_xlabel('log(Modulus Size)')
    ax2.set_ylabel('log(Memory)')
    ax2.set_title('Memory Complexity Analysis')
    ax2.legend()
    
    plt.savefig('complexity_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def crypto_analysis(df):
    print("\n=== ANALYSIS FOR CRYPTOGRAPHIC MODULI (256/512/1024/2048/4096 bits) ===")

    crypto_bits = [256, 512, 1024, 2048, 4096]
    crypto_df = df[df['N_bits'].isin(crypto_bits)]
    
    if crypto_df.empty:
        print("No data available for moduli 256/512/1024/2048/4096 bits")
        return
    
    # Calculate total attack time
    crypto_df['t_attack_s'] = crypto_df['t_recover_a_s'] + crypto_df['t_cma_s']
    
    for bits in crypto_bits:
        bit_data = crypto_df[crypto_df['N_bits'] == bits]
        if bit_data.empty:
            print(f"No data available for {bits} bits")
            continue
            
        print(f"\n--- Analysis for N={bits} bits ---")
        print(f"Success rate: {bit_data['success'].mean() * 100:.2f}%")
        print(f"Average 'recover a' time: {bit_data['t_recover_a_s'].mean():.4f} seconds")
        print(f"Average CMA time: {bit_data['t_cma_s'].mean():.4f} seconds")
        print(f"Average total attack time: {bit_data['t_attack_s'].mean():.4f} seconds")
        print(f"Average peak memory: {bit_data['peak_mem_kb'].mean() / 1024:.2f} MB")
        
        # Check if time and memory are within practical limits
        max_practical_time = 3600  # 1 hour
        max_practical_memory = 16 * 1024 * 1024  # 16 GB in KB
        
        avg_time = bit_data['t_attack_s'].mean()
        avg_memory = bit_data['peak_mem_kb'].mean()
        
        if avg_time > max_practical_time:
            print(f"⚠️  WARNING: Attack time ({avg_time:.2f} seconds) may not be practical")
        else:
            print(f"✅ Practical attack time: {avg_time:.2f} seconds")
            
        if avg_memory > max_practical_memory:
            print(f"⚠️  WARNING: Memory usage ({avg_memory/1024:.2f} MB) may not be practical")
        else:
            print(f"✅ Practical memory usage: {avg_memory/1024:.2f} MB")

def main():
    print("Scalability Analysis of Attacks on Cubic Pell Curve Cryptography Schemes")
    print("=" * 80)
    
    # Load data
    try:
        df = pd.read_csv('/home/dandi/smt8/papernew/eksperimen/cpc_attack_bench.csv')
    except FileNotFoundError:
        print("File 'cpc_attack_bench.csv' not found.")
        return
    
    success_rate = basic_analysis(df)
    runtime_stats = runtime_analysis(df)
    memory_stats = memory_analysis(df)
    complexity_results = complexity_analysis(df)
    crypto_analysis(df)
    
    # Create visualizations
    create_visualizations(df, complexity_results)
    
    print("\n=== CONCLUSION ===")
    print("The plots have been saved as:")
    print("- success_rate.png: Attack success rate")
    print("- attack_runtime.png: Attack operation runtime")
    print("- memory_usage.png: Memory usage")
    print("- complexity_analysis.png: Time and memory complexity analysis")
    
    time_complexity = complexity_results['time_complexity']
    mem_complexity = complexity_results['memory_complexity']
    
    print(f"\nBased on the log-log regression analysis:")
    print(f"- Attack time complexity: O(n^{time_complexity:.3f})")
    print(f"- Memory complexity: O(n^{mem_complexity:.3f})")

if __name__ == "__main__":
    main()
