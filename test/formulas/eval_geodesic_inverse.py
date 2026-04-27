import pandas
import matplotlib.pyplot as plt


def compute_metrics(data: pandas.DataFrame):
    for name in ['distance', 'azimuth', 'reverse_azimuth', 'reduced_length', 'geodesic_scale']:
        data[f'error_{name}'] = abs(data[f'boost_{name}'] - data[f'geographiclib_{name}'])
        data[f'relative_error_{name}'] = data[f'error_{name}'] / data[f'geographiclib_{name}']


def main():
    data_old = pandas.read_csv('results_old.csv', delimiter=',')
    compute_metrics(data_old)
    data_new = pandas.read_csv('results_new.csv', delimiter=',')
    compute_metrics(data_new)

    for name in ['distance', 'azimuth', 'reverse_azimuth', 'reduced_length', 'geodesic_scale']:
        _, axs = plt.subplots(2, 3, layout='compressed')
        axs[0, 0].loglog(data_old['geographiclib_distance'], data_old[f'error_{name}'], 'r.', alpha=0.1, markersize=2)
        axs[0, 0].set_ylabel(f'absolute error for {name}')
        axs[0, 1].loglog(data_new['geographiclib_distance'], data_new[f'error_{name}'], 'g.', alpha=0.1, markersize=2)
        axs[0, 1].sharey(axs[0, 0])
        axs[0, 2].loglog(data_old['geographiclib_distance'], data_old[f'error_{name}'], 'r.', alpha=0.1, markersize=2, label='before')
        axs[0, 2].loglog(data_new['geographiclib_distance'], data_new[f'error_{name}'], 'g.', alpha=0.1, markersize=2, label='after')
        axs[0, 2].sharey(axs[0, 0])
        legend = axs[0, 2].legend()
        for lh in legend.legend_handles:
            lh.set_alpha(1)

        axs[1, 0].loglog(data_old['geographiclib_distance'], data_old[f'relative_error_{name}'], 'r.', alpha=0.1, markersize=2)
        axs[1, 0].set_xlabel('distance')
        axs[1, 0].set_ylabel(f'relative error for {name}')
        axs[1, 1].loglog(data_new['geographiclib_distance'], data_new[f'relative_error_{name}'], 'g.', alpha=0.1, markersize=2)
        axs[1, 1].set_xlabel('distance')
        axs[1, 1].sharey(axs[1, 0])
        axs[1, 2].loglog(data_old['geographiclib_distance'], data_old[f'relative_error_{name}'], 'r.', alpha=0.1, markersize=2, label='before')
        axs[1, 2].loglog(data_new['geographiclib_distance'], data_new[f'relative_error_{name}'], 'g.', alpha=0.1, markersize=2, label='after')
        axs[1, 2].set_xlabel('distance')
        axs[1, 2].sharey(axs[1, 0])
        legend = axs[1, 2].legend()
        for lh in legend.legend_handles:
            lh.set_alpha(1)

        plt.show()


if __name__ == "__main__":
    main()
