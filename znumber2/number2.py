import numpy as np
import matplotlib.pyplot as plt

# ========== ПАРАМЕТРЫ ВАРИАНТА 18 ==========
f = 8.5e9            # частота, Гц
c = 3e8              # скорость света, м/с
lam = c / f          # длина волны, м
k = 2 * np.pi / lam  # волновое число
length_ratio = 0.5   # 2l / λ
l = (length_ratio / 2) * lam   # длина плеча, м

print(f"Частота: {f/1e9} ГГц")
print(f"Длина волны: {lam:.4f} м")
print(f"Длина плеча l: {l:.4f} м")

# ========== АНАЛИТИЧЕСКИЙ РАСЧЁТ ==========
theta = np.linspace(1e-6, np.pi - 1e-6, 1000)   # радианы
theta_deg = theta * 180 / np.pi
sin_theta = np.sin(theta)

E_theta = (np.cos(k * l * np.cos(theta)) - np.cos(k * l)) / sin_theta
F = np.abs(E_theta)
F = F / np.max(F)

integral = np.trapezoid(F**2 * sin_theta, theta)   # ∫₀^π F² sinθ dθ
D_max = 2 / integral                               # упрощённая формула
D_max_dB = 10 * np.log10(D_max)

print(f"\nМаксимальный КНД (аналитика):")
print(f"D_max = {D_max:.4f} (разы)")
print(f"D_max = {D_max_dB:.2f} дБ")

D_theta = D_max * F**2
D_theta_dB = 10 * np.log10(D_theta)

# ========== ЗАГРУЗКА ДАННЫХ ИЗ CST (rcs_results.txt) ==========
cst_file = "rcs_results.txt"
angles_deg = []
gains_dB = []

try:
    with open(cst_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                angle = float(parts[1])
                value_str = parts[2].replace(',', '.')
                value = float(value_str)
                angles_deg.append(angle)
                gains_dB.append(value)
    angles_deg = np.array(angles_deg)
    gains_dB = np.array(gains_dB)
    gains_lin = 10 ** (gains_dB / 10)
    print(f"\nЗагружено {len(angles_deg)} точек из {cst_file}")
except Exception as e:
    print(f"\nОшибка загрузки {cst_file}: {e}")
    angles_deg = None

# ========== ПОСТРОЕНИЕ ОБЩЕГО РИСУНКА С 4 ПОДГРАФИКАМИ ==========
fig = plt.figure(figsize=(12, 10))
fig.suptitle('Диаграмма направленности вибратора (вариант 18)', fontsize=14)

# 1) Декартова, разы
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(theta_deg, D_theta, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    ax1.plot(angles_deg, gains_lin, 'r--', linewidth=1.5, label='CST')
ax1.set_xlabel('θ, градусы')
ax1.set_ylabel('КНД, разы')
ax1.set_title('Декартова система (разы)')
ax1.grid(True)
ax1.legend()

# 2) Декартова, дБ
ax2 = fig.add_subplot(2, 2, 2)
ax2.plot(theta_deg, D_theta_dB, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    ax2.plot(angles_deg, gains_dB, 'r--', linewidth=1.5, label='CST')
ax2.set_xlabel('θ, градусы')
ax2.set_ylabel('КНД, дБ')
ax2.set_title('Декартова система (дБ)')
ax2.grid(True)
ax2.legend()

# 3) Полярная, разы
ax3 = fig.add_subplot(2, 2, 3, projection='polar')
ax3.plot(theta, D_theta, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    ax3.plot(np.deg2rad(angles_deg), gains_lin, 'r--', linewidth=1.5, label='CST')
ax3.set_title('Полярная система (разы)')
ax3.legend()
ax3.grid(True)

# 4) Полярная, дБ
ax4 = fig.add_subplot(2, 2, 4, projection='polar')
ax4.plot(theta, D_theta_dB, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    ax4.plot(np.deg2rad(angles_deg), gains_dB, 'r--', linewidth=1.5, label='CST')
ax4.set_title('Полярная система (дБ)')
ax4.legend()
ax4.grid(True)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('dipole_comparison_variant18.png', dpi=150)
plt.show()

# ========== СОХРАНЕНИЕ ТАБЛИЦЫ (аналитика) ==========
output_file = "dipole_results_variant18.txt"
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("theta_deg    D(разы)        D(дБ)\n")
    for deg, d_lin, d_db in zip(theta_deg, D_theta, D_theta_dB):
        f.write(f"{deg:8.4f}  {d_lin:10.6e}  {d_db:8.4f}\n")

print(f"\nАналитические данные сохранены в {output_file}")
print("\nПервые 20 строк результата:")
with open(output_file, 'r') as f:
    for i, line in enumerate(f):
        if i >= 20:
            break
        print(line.rstrip())