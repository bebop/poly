module.exports = {
  title: 'Poly',
  tagline: 'A Go library and command line utility for engineering organisms.',
  url: 'https://your-docusaurus-test-site.com',
  baseUrl: '/',
  favicon: 'img/favicon.ico',
  organizationName: 'TimothyStiles', // Usually your GitHub org/user name.
  projectName: 'Poly', // Usually your repo name.
  themeConfig: {
    navbar: {
      title: 'Poly',
      logo: {
        alt: 'My Site Logo',
        src: 'img/logo.svg',
      },
      links: [
        {
          to: 'docs/',
          activeBasePath: 'docs',
          label: 'Docs',
          position: 'left',
        },
        {to: 'blog', label: 'Blog', position: 'left'},
        {
          href: 'https://github.com/TimothyStiles/poly',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            // {
            //   label: 'Style Guide',
            //   to: 'docs/',
            // },
          ],
        },
        {
          title: 'Community',
          items: [
            // {
            //   label: 'Discord',
            //   href: '',
            // },
            {
              label: 'Twitter',
              href: 'https://twitter.com/timothystiles',
            },
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'Blog',
              to: 'blog',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/timothystiles/poly',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} Timothy Stiles`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          // It is recommended to set document id as docs home page (`docs/` path).
          homePageId: 'installation',
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/timothystiles/poly/edit/prime/',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/timothystiles/poly/edit/prime/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
