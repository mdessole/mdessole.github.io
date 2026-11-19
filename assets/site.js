// The content and navigation also work without JavaScript.
const publicationTools = document.querySelector('.publication-tools');
if (publicationTools) {
  const filters = [...publicationTools.querySelectorAll('[data-filter]')];
  const publications = [...document.querySelectorAll('.publication')];
  const count = document.querySelector('#publication-count');
  for (const button of filters) {
    button.addEventListener('click', () => {
      const type = button.dataset.filter;
      let visible = 0;
      for (const publication of publications) {
        publication.hidden = type !== 'all' && publication.dataset.type !== type;
        if (!publication.hidden) visible += 1;
      }
      for (const filter of filters) filter.setAttribute('aria-pressed', String(filter === button));
      count.textContent = `${visible} publication${visible === 1 ? '' : 's'}`;
    });
  }
  count.textContent = `${publications.length} publications`;
  publicationTools.hidden = false;
}
const year = document.querySelector('#year');
if (year) year.textContent = new Date().getFullYear();
