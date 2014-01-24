$(document).ready(function() {

    $(document).on('click',"div.collapsible div.unitview-header", function(){
        $(this).parent().find('.unitview-content').slideToggle(); 
        $(this).parent().find('.footer').slideToggle();
        $(this).parent().find('.unitview-counter').fadeToggle(); 
        $(this).parent().toggleClass(function() {
              if ($(this).is('.collapsed')) {
                return 'extended';
              } else {
                return 'collapsed';
              }
            });
    });


    


    $("div.collapsed .unitview-content").livequery(function(){
        $(this).hide()
    });

    $("div.collapsed .unitview-footer").livequery(function(){
        $(this).hide()
    });


//    DATATABLE HANDLING

    $('.datatable').each(function() {
            $(this).dataTable({
                "sDom": "<'row'<'span5'l><'span5'f>r>t<'row'<'span5'i><'span5'p>>",
                "sWrapper": "dataTables_wrapper form-inline",
                "sPaginationType": "bootstrap",
                "bPaginate": true,
                "bSort": true,
                "bLengthChange": true,
                "bFilter": true,
                "bInfo": true,
                "bAutoWidth": true ,
                "bDeferRender": true,
                "oLanguage": { "sSearch": "Filter:" }
            });
    });

    $('.datatable-not-sorted').livequery(function() {
        $(this).dataTable({
            "sDom": "<'row'<'span5'l><'span5'f>r>t<'row'<'span5'i><'span5'p>>",
            "sWrapper": "dataTables_wrapper form-inline",
            "sPaginationType": "bootstrap",
            "bPaginate": true,
            "bSort": true,
            "bLengthChange": true,
            "bFilter": true,
            "bInfo": true,
            "bAutoWidth": true ,
            "bDeferRender": true,
            "aaSorting": [],
            "oLanguage": { "sSearch": "Filter:" }

        });
    });

    $('.datatable-inverted').livequery(function() {
            $(this).dataTable({
                "sDom": "<'row'<'span5'l><'span5'f>r>t<'row'<'span5'i><'span5'p>>",
                "sWrapper": "dataTables_wrapper form-inline",
                "sPaginationType": "bootstrap",
                "bPaginate": true,
                "bSort": true,
                "bLengthChange": true,
                "bFilter": true,
                "bInfo": true,
                "bAutoWidth": true ,
                "bDeferRender": true,
                "aaSorting": [[ 0, "desc" ]],
                "oLanguage": { "sSearch": "Filter:" }
            });
    });
    
    



//    $('.editable').livequery(function () {
//        $(this).append('<span class="fg-button ui-state-default fg-button-icon-center ui-corner-all edit-icon"><span class="ui-icon ui-icon-pencil "></span></span>');
//        $(this).find('span.edit-icon').hide()
//     });
//
//    $('.editable').hover( function () {
//        $(this).find('span.edit-icon').show()},
//        function () {
//        $(this).find('span.edit-icon').hide();
//    });

//    $("div#view-feature-dialog").livequery(function(){
//        $(this).dialog({
//            autoOpen: false,
//            modal: true,
//            title: 'View feature on sequence',
//            width: 700,
//            //height: 500,
//            });
//    });

    // fix for datatable to work with bootstrap
    $.extend( $.fn.dataTableExt.oStdClasses, {
        "sSortAsc": "header headerSortDown",
        "sSortDesc": "header headerSortUp",
        "sSortable": "header"
    } );


    //scrollpsy
    $(".container").scrollspy();
    $('[data-spy="scroll"]').each(function () {
        var $spy = $(this).scrollspy('refresh')
    });




});
