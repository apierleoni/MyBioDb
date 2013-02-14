$(document).ready(function() {

    $("div.collapsible div.unitview-header").live('click', function(){
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



    $('.datatable').each(function() {
            $(this).dataTable({
                "sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
                "sWrapper": "dataTables_wrapper form-inline",
                "sPaginationType": "bootstrap",
                "bPaginate": true,
                "bSort": true,
                "bLengthChange": true,
                "bFilter": true,
                "bInfo": true,
                "bAutoWidth": true ,
                "bDeferRender": true,
            });
    });

    $('.datatable-not-sorted').livequery(function() {
        $(this).dataTable({
            "sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
            "sWrapper": "dataTables_wrapper form-inline",
            "sPaginationType": "bootstrap",
            "bPaginate": true,
            "bSort": true,
            "bLengthChange": true,
            "bFilter": true,
            "bInfo": true,
            "bAutoWidth": true ,
            "sPaginationType": "full_numbers",
            "bDeferRender": true,
            "aaSorting": []
        });
    });

    $('.datatable-inverted').livequery(function() {
            $(this).dataTable({
                "sDom": "<'row'<'span6'l><'span6'f>r>t<'row'<'span6'i><'span6'p>>",
                "sWrapper": "dataTables_wrapper form-inline",
                "sPaginationType": "bootstrap",
                "bPaginate": true,
                "bSort": true,
                "bLengthChange": true,
                "bFilter": true,
                "bInfo": true,
                "bAutoWidth": true ,
                "bDeferRender": true,
                "aaSorting": [[ 0, "desc" ]]
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

});
