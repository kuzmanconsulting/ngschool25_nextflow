document.addEventListener("DOMContentLoaded", function() {
    hljs.registerLanguage("nextflow", function(hljs) {
        return {
            keywords: {
                keyword: 'process script input output params channel publishDir tuple'
            },
            contains: [
                hljs.COMMENT('//', '$'),
                hljs.COMMENT('/\\*', '\\*/'),
                {
                    className: 'string',
                    begin: '"', end: '"'
                },
                {
                    className: 'number',
                    begin: '\\b\\d+(\\.\\d+)?'
                },
                {
                    className: 'built_in',
                    begin: '\\b(Channel|file|val)\\b'
                }
            ]
        };
    });

    hljs.highlightAll(); // Reapply highlighting
});