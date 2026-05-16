import numpy
import pandas
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['figure.figsize'] = (19.20, 10.80)


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
        fig, axs = plt.subplots(2, 3, layout='compressed')
        axs[0, 0].loglog(data_old['geographiclib_distance'], data_old[f'error_{name}'], 'r.', alpha=0.1, markersize=2)
        axs[0, 0].set_ylabel(f'absolute error for {name}')
        axs[0, 1].loglog(data_new['geographiclib_distance'], data_new[f'error_{name}'], 'g.', alpha=0.1, markersize=2)
        axs[0, 1].sharey(axs[0, 0])
        axs[0, 2].loglog(data_old['geographiclib_distance'], data_old[f'error_{name}'], 'r.', alpha=0.1, markersize=2, label='before')
        axs[0, 2].loglog(data_new['geographiclib_distance'], data_new[f'error_{name}'], 'g.', alpha=0.1, markersize=2, label='after')
        axs[0, 2].sharey(axs[0, 0])

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
        fig.savefig(f'distsorted_{name}.png')

        fig, axs = plt.subplots(2, 2, layout='compressed')

        perm = numpy.argsort(data_new[f'error_{name}'])
        axs[0, 0].semilogy(data_old[f'error_{name}'].values[perm], 'r.', alpha=0.5, markersize=2, label='before')
        axs[0, 0].semilogy(data_new[f'error_{name}'].values[perm], 'g.', alpha=0.5, markersize=2, label='after')
        axs[0, 0].set_ylabel(f'absolute error for {name}')
        perm = numpy.argsort(data_old[f'error_{name}'])
        axs[0, 1].semilogy(data_new[f'error_{name}'].values[perm], 'g.', alpha=0.5, markersize=2, label='after')
        axs[0, 1].semilogy(data_old[f'error_{name}'].values[perm], 'r.', alpha=0.5, markersize=2, label='before')

        perm = numpy.argsort(data_new[f'relative_error_{name}'])
        axs[1, 0].semilogy(data_old[f'relative_error_{name}'].values[perm], 'r.', alpha=0.5, markersize=2, label='before')
        axs[1, 0].semilogy(data_new[f'relative_error_{name}'].values[perm], 'g.', alpha=0.5, markersize=2, label='after')
        axs[1, 0].set_ylabel(f'relative error for {name}')
        perm = numpy.argsort(data_old[f'relative_error_{name}'])
        axs[1, 1].semilogy(data_new[f'relative_error_{name}'].values[perm], 'g.', alpha=0.5, markersize=2, label='after')
        axs[1, 1].semilogy(data_old[f'relative_error_{name}'].values[perm], 'r.', alpha=0.5, markersize=2, label='before')
        for lh in axs[1, 1].legend().legend_handles:
            lh.set_alpha(1)
        fig.savefig(f'errsorted_{name}.png')


if __name__ == "__main__":
    main()
