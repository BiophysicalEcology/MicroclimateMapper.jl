import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}
const navTemp = {
  nav: [
    { text: 'Home', link: '/' },
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Background',
          items: [
            { text: 'Introduction', link: '/manual/introduction' },
            { text: 'One keyword per data set', link: '/manual/data_sources' },
            { text: 'Built on Rasters and DimensionalData', link: '/manual/spatial_stack' },
          ]
        },
        { text: 'Data',
          items: [
            { text: 'Forcing data sets', link: '/manual/forcing_data' },
            { text: 'Terrain and atmosphere', link: '/manual/terrain' },
            { text: 'Surfaces and soils', link: '/manual/surfaces_soils' },
          ]
        },
        { text: 'Running',
          items: [
            { text: 'Outputs', link: '/manual/outputs' },
            { text: 'Performance', link: '/manual/performance' },
            { text: 'Lateral flows', link: '/manual/lateral_flows' },
          ]
        },
        { text: 'Reference',
          items: [
            { text: 'For NicheMapR users', link: '/manual/nichemapr' },
            { text: 'References', link: '/manual/references' },
          ]
        },
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'Points', link: '/tutorials/points' },
        { text: 'Swapping data sets', link: '/tutorials/swap_data' },
        { text: 'Maps', link: '/tutorials/maps' },
        { text: 'Terrain: Saba', link: '/tutorials/saba' },
        { text: 'Canopy on Saba', link: '/tutorials/canopy' },
        { text: 'Soils and land cover', link: '/tutorials/soils_landcover' },
      ]
    },
    { text: 'Interactive', link: '/interactive' },
    { text: 'API', link: '/api' }
  ],
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "Microclimates over points and grids, from any spatial data set, with Microclimate.jl",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  ignoreDeadLinks: false,
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },

  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
    { text: 'Get Started', link: '/get_started' },
    { text: 'Manual',
      items: [
        { text: 'Background', collapsed: false,
          items: [
            { text: 'Introduction', link: '/manual/introduction' },
            { text: 'One keyword per data set', link: '/manual/data_sources' },
            { text: 'Built on Rasters and DimensionalData', link: '/manual/spatial_stack' },
          ]
        },
        { text: 'Data', collapsed: false,
          items: [
            { text: 'Forcing data sets', link: '/manual/forcing_data' },
            { text: 'Terrain and atmosphere', link: '/manual/terrain' },
            { text: 'Surfaces and soils', link: '/manual/surfaces_soils' },
          ]
        },
        { text: 'Running', collapsed: false,
          items: [
            { text: 'Outputs', link: '/manual/outputs' },
            { text: 'Performance', link: '/manual/performance' },
            { text: 'Lateral flows', link: '/manual/lateral_flows' },
          ]
        },
        { text: 'Reference', collapsed: false,
          items: [
            { text: 'For NicheMapR users', link: '/manual/nichemapr' },
            { text: 'References', link: '/manual/references' },
          ]
        },
      ]
    },
    { text: 'Tutorials',
      items: [
        { text: 'Points', link: '/tutorials/points' },
        { text: 'Swapping data sets', link: '/tutorials/swap_data' },
        { text: 'Maps', link: '/tutorials/maps' },
        { text: 'Terrain: Saba', link: '/tutorials/saba' },
        { text: 'Canopy on Saba', link: '/tutorials/canopy' },
        { text: 'Soils and land cover', link: '/tutorials/soils_landcover' },
      ]
    },
    { text: 'Interactive', link: '/interactive' },
    { text: 'API', link: '/api' }
    ],
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/" target="_blank"><strong>DocumenterVitepress.jl</strong></a> <br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}. Released under the MIT License.`
    }
  }
})
