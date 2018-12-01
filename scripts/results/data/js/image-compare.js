(function($){
    $.fn.extend({
        //plugin name - qcrosscompare
        qcrosscompare: function(options) {

            return this.each(function() {

	            var i = $(this);
				var img_layer_tl = i.children('img:eq(0)').attr('src');
				var img_layer_tr = i.children('img:eq(1)').attr('src');
				var img_layer_bl = i.children('img:eq(2)').attr('src');
				var img_layer_br = i.children('img:eq(3)').attr('src');
				var img_cap_tl = i.children('img:eq(0)').attr('alt');
				var img_cap_tr = i.children('img:eq(1)').attr('alt');
				var img_cap_bl = i.children('img:eq(2)').attr('alt');
				var img_cap_br = i.children('img:eq(3)').attr('alt');

				var width = i.children('img:eq(0)').width();
				var height = i.children('img:eq(0)').height();

				i.children('img').hide();

				i.css({'overflow': 'hidden', 'position': 'relative'});
				i.append('<div class="ba-layer_tl"></div>');
				i.append('<div class="ba-layer_tr"></div>');
				i.append('<div class="ba-layer_bl"></div>');
				i.append('<div class="ba-layer_br"></div>');

				i.append('<div class="ba-caption_tl">' + img_cap_tl + '</div>');
				i.append('<div class="ba-caption_tr">' + img_cap_tr + '</div>');
				i.append('<div class="ba-caption_bl">' + img_cap_bl + '</div>');
				i.append('<div class="ba-caption_br">' + img_cap_br + '</div>');

				i.children('.ba-layer_tl, .ba-layer_tr, .ba-layer_bl, .ba-layer_br').width(width);
				i.children('.ba-layer_tl, .ba-layer_tr, .ba-layer_bl, .ba-layer_br').height(height);

				i.children('.ba-layer_tl').css('backgroundImage','url(' + img_layer_tl + ')');
				i.children('.ba-layer_tr').css('backgroundImage','url(' + img_layer_tr + ')');
				i.children('.ba-layer_bl').css('backgroundImage','url(' + img_layer_bl + ')');
				i.children('.ba-layer_br').css('backgroundImage','url(' + img_layer_br + ')');
				i.children('.ba-layer_tl').css('backgroundSize', width + 'px ' + height + 'px');
				i.children('.ba-layer_tr').css('backgroundSize', width + 'px ' + height + 'px');
				i.children('.ba-layer_bl').css('backgroundSize', width + 'px ' + height + 'px');
				i.children('.ba-layer_br').css('backgroundSize', width + 'px ' + height + 'px');

				i.children('.ba-layer_tl').width(width * 0.5);
				i.children('.ba-layer_tl').height(height * 0.5);
				i.children('.ba-layer_tr').height(height * 0.5);
				i.children('.ba-layer_bl').width(width * 0.5);

				i.children('.ba-caption_tl').show();
				i.children('.ba-caption_tr').show();
				i.children('.ba-caption_bl').show();
				i.children('.ba-caption_br').show();

				i.children('.ba-caption_tl').css({ bottom: height * 0.5, right: width * 0.5 });
				i.children('.ba-caption_tr').css({ bottom: height * 0.5, left:  width * 0.5 });
				i.children('.ba-caption_bl').css({ top:    height * 0.5, right: width * 0.5 });
				i.children('.ba-caption_br').css({ top:    height * 0.5, left:  width * 0.5 });

            }).mousemove(function (e) {

				var o = options;
				var i = $(this);

				right_border_width = parseInt(i.children('.ba-layer_tl').css("borderRightWidth"), 10);
				bottom_border_width = parseInt(i.children('.ba-layer_tl').css("borderBottomWidth"), 10);
				pos_imgX = this.getBoundingClientRect().left;
				pos_imgY = this.getBoundingClientRect().top;
				pos_mouseX = e.clientX - right_border_width * 0.5;
				pos_mouseY = e.clientY - bottom_border_width * 0.5;
				new_width  = pos_mouseX - pos_imgX;
				new_height = pos_mouseY - pos_imgY;
				img_width  = i.width();
				img_height = i.height();

				i.children('.ba-layer_tl').width(new_width);
				i.children('.ba-layer_tl').height(new_height);
				i.children('.ba-layer_tr').height(new_height);
				i.children('.ba-layer_bl').width(new_width);

				i.children('.ba-caption_tl').css({ bottom: img_height - new_height, right: img_width - new_width });
				i.children('.ba-caption_tr').css({ bottom: img_height - new_height, left:  new_width });
				i.children('.ba-caption_bl').css({ top:    new_height,              right: img_width - new_width });
				i.children('.ba-caption_br').css({ top:    new_height,              left:  new_width });
			});
        }
    });
	})(jQuery);

$(function () {
	$('.cross_compare').qcrosscompare();
});