pip install numpy matplotlib scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

# ========== ПАРАМЕТРЫ ВАРИАНТА 18 ==========
f = 8.5e9          # частота, Гц
c = 3e8            # скорость света, м/с
lam = c / f        # длина волны, м
k = 2 * np.pi / lam
length_ratio = 0.5          # 2l/λ
l = (length_ratio / 2) * lam   # длина плеча, м

print(f"Частота: {f/1e9} ГГц")
print(f"Длина волны: {lam:.4f} м")
print(f"Длина плеча l: {l:.4f} м")

# ========== АНАЛИТИЧЕСКИЙ РАСЧЁТ ==========
theta = np.linspace(1e-6, np.pi - 1e-6, 1000)   # радианы
theta_deg = theta * 180 / np.pi                  # градусы
sin_theta = np.sin(theta)

# Характеристика направленности по полю (без постоянных множителей)
E_theta = (np.cos(k * l * np.cos(theta)) - np.cos(k * l)) / sin_theta
F = np.abs(E_theta)
F /= np.max(F)                                   # нормировка

# Вычисление Dmax
integrand = F**2 * sin_theta
integral = simps(integrand, theta)               # ∫₀^π F² sinθ dθ
D_max = 2 / integral
D_max_dB = 10 * np.log10(D_max)

print(f"\nМаксимальный КНД (аналитика):")
print(f"D_max = {D_max:.4f} (разы)")
print(f"D_max = {D_max_dB:.2f} дБ")

# КНД в зависимости от угла
D_theta = D_max * F**2
D_theta_dB = 10 * np.log10(D_theta)

# ========== ЗАГРУЗКА ДАННЫХ ИЗ CST ==========
cst_file = "rcs_results.txt"
try:
    # Читаем файл, заменяем запятые на точки, извлекаем углы и дБ
    angles_deg = []
    gains_dB = []
    with open(cst_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                # Столбцы: номер, угол(град), значение(дБ)
                angle = float(parts[1])
                # Заменяем запятую на точку и преобразуем в float
                value_str = parts[2].replace(',', '.')
                value = float(value_str)
                angles_deg.append(angle)
                gains_dB.append(value)
    angles_deg = np.array(angles_deg)
    gains_dB = np.array(gains_dB)
    gains_lin = 10 ** (gains_dB / 10)
    print(f"\nЗагружены данные из {cst_file}, точек: {len(angles_deg)}")
    print(f"Диапазон углов: {angles_deg[0]}° ... {angles_deg[-1]}°")
    print(f"Диапазон КНД (дБ): {gains_dB.min():.2f} ... {gains_dB.max():.2f}")
except Exception as e:
    print(f"\nОшибка загрузки {cst_file}: {e}")
    print("Будут построены только аналитические кривые.")
    angles_deg = None
    gains_dB = None
    gains_lin = None

# ========== ПОСТРОЕНИЕ ГРАФИКОВ ==========
# 1. Декартова система, разы
plt.figure(figsize=(8, 5))
plt.plot(theta_deg, D_theta, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    plt.plot(angles_deg, gains_lin, 'r--', linewidth=1.5, label='CST (разы)')
plt.xlabel('Угол θ, градусы')
plt.ylabel('КНД D(θ), разы')
plt.title('Диаграмма направленности (декартова, разы)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig('D_cartesian_lin.png', dpi=150)
plt.show()

# 2. Декартова система, дБ
plt.figure(figsize=(8, 5))
plt.plot(theta_deg, D_theta_dB, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    plt.plot(angles_deg, gains_dB, 'r--', linewidth=1.5, label='CST (дБ)')
plt.xlabel('Угол θ, градусы')
plt.ylabel('КНД D(θ), дБ')
plt.title('Диаграмма направленности (декартова, дБ)')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig('D_cartesian_dB.png', dpi=150)
plt.show()

# 3. Полярная система, разы
plt.figure(figsize=(7, 7))
ax = plt.subplot(111, projection='polar')
ax.plot(theta, D_theta, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    angles_rad = np.deg2rad(angles_deg)
    ax.plot(angles_rad, gains_lin, 'r--', linewidth=1.5, label='CST (разы)')
ax.set_title('Диаграмма направленности (полярная, разы)')
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.savefig('D_polar_lin.png', dpi=150)
plt.show()

# 4. Полярная система, дБ
plt.figure(figsize=(7, 7))
ax = plt.subplot(111, projection='polar')
ax.plot(theta, D_theta_dB, 'b-', linewidth=2, label='Аналитика')
if angles_deg is not None:
    ax.plot(angles_rad, gains_dB, 'r--', linewidth=1.5, label='CST (дБ)')
ax.set_title('Диаграмма направленности (полярная, дБ)')
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.savefig('D_polar_dB.png', dpi=150)
plt.show()

# ========== СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ==========
output_file = "dipole_results_variant18.txt"
with open(output_file, 'w', encoding='utf-8') as f:
    f.write("theta(deg)  D(разы)    D(дБ)\n")
    for deg, d_lin, d_db in zip(theta_deg, D_theta, D_theta_dB):
        f.write(f"{deg:8.4f}  {d_lin:10.6e}  {d_db:8.4f}\n")

print(f"\nРезультаты (аналитика) сохранены в файл {output_file}")
print("\nПервые 20 строк аналитических данных:")
with open(output_file, 'r', encoding='utf-8') as f:
    for i, line in enumerate(f):
        if i >= 20:
            break
        print(line.rstrip())