import numpy as np
import matplotlib.pyplot as plt
import toml
import json
from scipy.special import spherical_jn, spherical_yn

c = 3e8  # скорость света (м/с)


# ---------- 1. Класс загрузки данных ----------
class VariantLoader:
    def __init__(self, filename, variant):
        self.filename = filename
        self.variant = variant

    def load(self):
        """
        Читает файл .toml и возвращает параметры нужного варианта
        """
        data = toml.load(self.filename)["data"]
        for item in data:
            if item["variant"] == self.variant:
                return {
                    "D": float(item["D"]),      # диаметр сферы
                    "fmin": float(item["fmin"]),# минимальная частота
                    "fmax": float(item["fmax"]) # максимальная частота
                }
        raise ValueError("Вариант не найден")


# ---------- 2. Класс расчета ЭПР ----------
class RCSCalculator:
    def __init__(self, r):
        self.r = r  # радиус сферы

    def hankel(self, n, x):
        """
        Сферическая функция Ханкеля:
        h_n = j_n + i*y_n
        """
        return spherical_jn(n, x) + 1j * spherical_yn(n, x)

    def calc_sigma(self, f, n_max=20):
        """
        Расчет ЭПР для одной частоты
        """
        lam = c / f                  # длина волны
        k = 2 * np.pi / lam          # волновое число
        kr = k * self.r              # аргумент функций

        s = 0  # сумма ряда

        for n in range(1, n_max + 1):
            # функции Бесселя
            jn = spherical_jn(n, kr)
            jn_1 = spherical_jn(n - 1, kr)

            # функции Ханкеля
            hn = self.hankel(n, kr)
            hn_1 = self.hankel(n - 1, kr)

            # коэффициенты
            a_n = jn / hn
            b_n = (kr * jn_1 - n * jn) / (kr * hn_1 - n * hn)

            # суммирование
            s += ((-1)**n) * (n + 0.5) * (b_n - a_n)

        # итоговая формула ЭПР
        sigma = (lam**2 / np.pi) * abs(s)**2
        return sigma.real


# ---------- 3. Класс сохранения ----------
class ResultWriter:
    def __init__(self, filename):
        self.filename = filename

    def write(self, freq, lam, sigma):
        """
        Сохраняет данные в JSON (формат 3)
        """
        data = {
            "freq": freq.tolist(),
            "lambda": lam.tolist(),
            "rcs": sigma.tolist()
        }

        with open(self.filename, "w") as f:
            json.dump(data, f, indent=4)


# ---------- 4. Класс построения графика ----------
class Plotter:
    @staticmethod
    def plot(freq, sigma):
        """
        Строит график зависимости ЭПР от частоты
        """
        plt.figure(figsize=(8, 5))

        # сам график
        plt.plot(freq / 1e9, sigma)

        # подписи
        plt.xlabel("Частота, ГГц")
        plt.ylabel("ЭПР, м²")
        plt.title("Зависимость ЭПР идеально проводящей сферы от частоты (вариант 18)")

        plt.grid(True)
        plt.tight_layout()
        plt.savefig("plotvar18.png")
        plt.show()


# ---------- MAIN ----------
def main():
    # 1. Загружаем параметры варианта 18
    loader = VariantLoader("task_rcs_01.toml", 18)
    params = loader.load()

    r = params["D"] / 2  # радиус сферы
    fmin = params["fmin"]
    fmax = params["fmax"]

    # 2. Создаем массив частот
    freq = np.linspace(fmin, fmax, 400)
    lam = c / freq

    # 3. Расчет ЭПР
    calc = RCSCalculator(r)
    sigma = np.array([calc.calc_sigma(f) for f in freq])

    # 4. Сохранение результата
    writer = ResultWriter("rcs_18.json")
    writer.write(freq, lam, sigma)

    # 5. Вывод первых 20 строк (для отчета)
    print("\nПервые 20 строк результата:")
    print("freq (Гц)        lambda (м)        rcs (м^2)")
    for i in range(20):
        print(f"{freq[i]:.6e}   {lam[i]:.6e}   {sigma[i]:.6e}")

    # 6. Построение графика
    Plotter.plot(freq, sigma)


if __name__ == "__main__":
    main()